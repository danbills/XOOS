#include "segmentation/genomic-segments.h"

#include <algorithm>
#include <armadillo>
#include <cmath>
#include <set>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <csv.hpp>

#include <xoos/log/logging.h>

#include "likelihood/likelihood-flags.h"
#include "segmentation/genomic-segment-column-names.h"
#include "segmentation/interval-trees.h"
#include "utility/utility-functions.h"

namespace xoos::cnc::segmentation {

/**
 * @brief constructor for GenomicSegment
 * @param segment a Segment object to be converted to GenomicSegmnt
 * @param obvs the corresponding observations for the Segmnet object
 */
GenomicSegment::GenomicSegment(const Segment& segment, const Observations& obvs) {
  contig = obvs.contigs[segment.start];
  start = obvs.starts[segment.start];
  // segment.end is 0-based exclusive, so we have to substract one to access the last position
  end = obvs.ends[segment.end - 1];
  arr_start = segment.start;
  arr_end = segment.end;
}

/**
 * @brief populate the .mean_logr and the .num_obs fields
 * @param obvs corresponding observations
 */
void GenomicSegment::PopulateOptionalFields(const Observations& obvs) {
  arma::vec obvs_vec = obvs.obvs(arma::span(arr_start.value(), arr_end.value() - 1));
  mean_logr = arma::mean(obvs_vec);
  num_obs = obvs_vec.n_elem;
}

void GenomicSegment::SetField(const std::string& field, const std::string& val) {
  if (field == kGenomicSegColId) {
    id = val;
  } else if (field == kGenomicSegColContig) {
    contig = val;
  } else if (field == kGenomicSegColStart) {
    // the file format specifies 1-based, but we convert it to 0-based closed here for internal consistency
    start = StringToNonNegativeUIntOrThrow(val) - 1;
  } else if (field == kGenomicSegColEnd) {
    // the file format specifies 1-based, but we convert it to 0-based open here for internal consistency
    end = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColArrStart) {
    arr_start = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColArrEnd) {
    arr_end = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColNumLogr) {
    num_obs = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColMeanLogr) {
    mean_logr = std::stod(val);
  } else if (field == kGenomicSegColNumSnp) {
    num_snps = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColMBaf) {
    if (val == "nan") {
      mbaf = std::nullopt;
    } else {
      mbaf = std::stod(val);
    }
  } else if (field == kGenomicSegColMeanDh) {
    mean_dh = std::stod(val);
  } else if (field == kGenomicSegColAvgMeanMapq) {
    if (val == "NA" || val == "nan") {
      avg_mean_mapq = std::nullopt;
    } else {
      avg_mean_mapq = std::stod(val);
    }
  } else if (field == kGenomicSegColPloidy) {
    ploidy = std::stod(val);
  } else if (field == kGenomicSegColPurity) {
    purity = std::stod(val);
  } else if (field == kGenomicSegColTotalCopyNumber) {
    total_copy_number = StringToNonNegativeUIntOrThrow(val);
  } else if (field == kGenomicSegColExpectedTotalCopyNumber) {
    expected_total_copy_number = std::stod(val);
  } else if (field == kGenomicSegColMinorCopyNumber) {
    if (val == "nan") {
      minor_copy_number = std::nullopt;
    } else {
      minor_copy_number = StringToNonNegativeUIntOrThrow(val);
    }
  } else if (field == kGenomicSegColMajorCopyNumber) {
    if (val == "nan") {
      major_copy_number = std::nullopt;
    } else {
      major_copy_number = StringToNonNegativeUIntOrThrow(val);
    }
  } else if (field == kGenomicSegColLogrLikelihood) {
    if (val == "nan") {
      logr_likelihood = std::nullopt;
    } else {
      logr_likelihood = std::stod(val);
    }
  } else if (field == kGenomicSegColBafLikelihood) {
    if (val == "nan") {
      baf_likelihood = std::nullopt;
    } else {
      baf_likelihood = std::stod(val);
    }
  } else if (field == kGenomicSegColJointLikelihood) {
    if (val == "nan") {
      joint_likelihood = std::nullopt;
    } else {
      joint_likelihood = std::stod(val);
    }
  } else if (field == kGenomicSegColSex) {
    sex = ParseSexOption(val);
  } else if (field == kGenomicSegColFilter) {
    flags = StringToFlags(val);
  } else if (field == kGenomicSegColInAllosome) {
    if (val == "true") {
      in_allosome = true;
    } else if (val == "false") {
      in_allosome = false;
    } else {
      throw std::runtime_error("invalid value for in_allosome");
    }
  } else if (field == kGenomicSegColInPseudoAutosomalRegion) {
    if (val == "true") {
      in_pseudo_autosomal_region = true;
    } else if (val == "false") {
      in_pseudo_autosomal_region = false;
    } else {
      throw std::runtime_error("invalid value for in_pseudo_autosomal_region");
    }
  } else {
    Logging::Error("GenomicSegCol {} not supported", field);
    throw std::runtime_error("unsupported column");
  }
}

void GenomicSegment::SetFieldOrThrowError(const std::string& field, const std::string& val) {
  try {
    SetField(field, val);
  } catch (std::invalid_argument& e) {
    Logging::Error("Invalid value '{}' for field '{}': {}", val, field, e.what());
    throw std::runtime_error("invalid value for field");
  } catch (std::out_of_range& e) {
    Logging::Error("Out of range value '{}' for field '{}': {}", val, field, e.what());
    throw std::runtime_error("out of range value for field");
  } catch (std::runtime_error& e) {
    Logging::Error("Runtime error for field '{}': {}", field, e.what());
    throw;
  }
}

std::string GenomicSegment::GetFieldAsString(const std::string& field) const {
  if (field == kGenomicSegColId) {
    return id.value();
  } else if (field == kGenomicSegColContig) {
    return contig;
  } else if (field == kGenomicSegColStart) {
    return fmt::format("{}", start);
  } else if (field == kGenomicSegColEnd) {
    return fmt::format("{}", end);
  } else if (field == kGenomicSegColArrStart) {
    return fmt::format("{}", arr_start.value());
  } else if (field == kGenomicSegColArrEnd) {
    return fmt::format("{}", arr_end.value());
  } else if (field == kGenomicSegColNumLogr) {
    return fmt::format("{}", num_obs.value());
  } else if (field == kGenomicSegColMeanLogr) {
    return fmt::format("{:.2f}", mean_logr.value());
  } else if (field == kGenomicSegColNumSnp) {
    return fmt::format("{}", num_snps.value());
  } else if (field == kGenomicSegColMBaf) {
    if (!mbaf.has_value()) {
      return "nan";
    } else {
      return fmt::format("{:.2f}", mbaf.value());
    }
  } else if (field == kGenomicSegColMeanDh) {
    return fmt::format("{:.2f}", mean_dh.value());
  } else if (field == kGenomicSegColAvgMeanMapq) {
    if (!avg_mean_mapq.has_value()) {
      return "NA";
    }
    return fmt::format("{:.2f}", avg_mean_mapq.value());
  } else if (field == kGenomicSegColPloidy) {
    return fmt::format("{:.2f}", ploidy.value());
  } else if (field == kGenomicSegColPurity) {
    return fmt::format("{:.2f}", purity.value());
  } else if (field == kGenomicSegColTotalCopyNumber) {
    return fmt::format("{}", total_copy_number.value());
  } else if (field == kGenomicSegColExpectedTotalCopyNumber) {
    return fmt::format("{:.2f}", expected_total_copy_number.value());
  } else if (field == kGenomicSegColMinorCopyNumber) {
    if (!minor_copy_number.has_value() && !major_copy_number.has_value()) {
      return "nan";
    } else {
      return fmt::format("{}", minor_copy_number.value());
    }
  } else if (field == kGenomicSegColMajorCopyNumber) {
    if (!minor_copy_number.has_value() && !major_copy_number.has_value()) {
      return "nan";
    } else {
      return fmt::format("{}", major_copy_number.value());
    }
  } else if (field == kGenomicSegColLogrLikelihood) {
    if (!logr_likelihood.has_value()) {
      return "nan";
    }
    return fmt::format("{:.2f}", logr_likelihood.value());
  } else if (field == kGenomicSegColBafLikelihood) {
    if (!baf_likelihood.has_value()) {
      return "nan";
    }
    return fmt::format("{:.2f}", baf_likelihood.value());
  } else if (field == kGenomicSegColJointLikelihood) {
    if (!joint_likelihood.has_value()) {
      return "nan";
    }
    return fmt::format("{:.2f}", joint_likelihood.value());
  } else if (field == kGenomicSegColSex) {
    return SexToStr(sex.value());
  } else if (field == kGenomicSegColFilter) {
    return FlagsToString(flags.value());
  } else if (field == kGenomicSegColInAllosome) {
    return in_allosome.value() ? "true" : "false";
  } else if (field == kGenomicSegColInPseudoAutosomalRegion) {
    return in_pseudo_autosomal_region.value() ? "true" : "false";
  } else {
    Logging::Error("{} is an unsupported SEG column", field);
    throw std::runtime_error("invalid column");
  }
}

std::vector<std::string> GenomicSegment::GetPopulatedOptionalCols() const {
  std::vector<std::string> ret;
  if (arr_start.has_value()) {
    ret.emplace_back(kGenomicSegColArrStart);
  }
  if (arr_end.has_value()) {
    ret.emplace_back(kGenomicSegColArrEnd);
  }
  if (id.has_value()) {
    ret.emplace_back(kGenomicSegColId);
  }
  if (num_obs.has_value()) {
    ret.emplace_back(kGenomicSegColNumLogr);
  }
  if (num_snps.has_value()) {
    ret.emplace_back(kGenomicSegColNumSnp);
  }
  if (avg_mean_mapq.has_value()) {
    ret.emplace_back(kGenomicSegColAvgMeanMapq);
  }
  if (ploidy.has_value()) {
    ret.emplace_back(kGenomicSegColPloidy);
  }
  if (purity.has_value()) {
    ret.emplace_back(kGenomicSegColPurity);
  }
  if (total_copy_number.has_value()) {
    ret.emplace_back(kGenomicSegColTotalCopyNumber);
  }
  if (expected_total_copy_number.has_value()) {
    ret.emplace_back(kGenomicSegColExpectedTotalCopyNumber);
  }
  if (minor_copy_number.has_value()) {
    ret.emplace_back(kGenomicSegColMinorCopyNumber);
  }
  if (major_copy_number.has_value()) {
    ret.emplace_back(kGenomicSegColMajorCopyNumber);
  }
  if (logr_likelihood.has_value()) {
    ret.emplace_back(kGenomicSegColLogrLikelihood);
  }
  if (baf_likelihood.has_value()) {
    ret.emplace_back(kGenomicSegColBafLikelihood);
  }
  if (joint_likelihood.has_value()) {
    ret.emplace_back(kGenomicSegColJointLikelihood);
  }
  if (sex.has_value()) {
    ret.emplace_back(kGenomicSegColSex);
  }
  if (flags.has_value()) {
    ret.emplace_back(kGenomicSegColFilter);
  }
  if (mbaf.has_value()) {
    ret.emplace_back(kGenomicSegColMBaf);
  }
  if (in_allosome.has_value()) {
    ret.emplace_back(kGenomicSegColInAllosome);
  }
  if (in_pseudo_autosomal_region.has_value()) {
    ret.emplace_back(kGenomicSegColInPseudoAutosomalRegion);
  }
  if (mean_dh.has_value()) {
    ret.emplace_back(kGenomicSegColMeanDh);
  }
  if (mean_logr.has_value()) {
    ret.emplace_back(kGenomicSegColMeanLogr);
  }
  return ret;
}

/**
 * @brief converts a Segment to a GenomicSegment (i.e. translates it into genomic coordinates). Requires original
 * observations to complete the conversion.
 * TODO: require a `.genome` file so that we can extend segments to cover whole genome
 * @param obvs original observations from which `segments` was produced
 * @param segments segments generated from obvs
 * @return std::vector<GenomicSegment>, which are `segments` translated to genomic coordinates
 */
std::vector<GenomicSegment> ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                             const std::vector<Segment>& segments) {
  std::vector<GenomicSegment> genomic_segments;
  genomic_segments.reserve(segments.size());
  for (const auto& segment : segments) {
    genomic_segments.emplace_back(segment, obvs);
  }
  return genomic_segments;
}

/**
 * @brief same as above function, but appends results to a given std::vector
 * TODO: require a `.genome` file so that we can extend segments to cover whole genome
 * @param obvs original observations from which `segments` was produced
 * @param segments segments generated from obvs
 * @return std::vector<GenomicSegment>, which are `segments` translated to genomic coordinates
 */
std::vector<GenomicSegment>& ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                              const std::vector<Segment>& segments,
                                                              std::vector<GenomicSegment>& genomic_segments) {
  for (const auto& segment : segments) {
    genomic_segments.emplace_back(segment, obvs);
  }
  return genomic_segments;
}

/**
 * @brief adjusts the start and ends fields in each GenomicSegment in a list by setting the end to the midpoint of the
 * current GenomicSegment and the next, and setting the start to the midpoint of the current GenomicSegment and the
 * previous
 * @param genomic_segments std::vector of GenomicSegments
 * @return modifies the given GenomicSegments object
 */
std::vector<GenomicSegment>& AdjustGenomicSegmentBreakpoints(std::vector<GenomicSegment>& genomic_segments) {
  for (size_t i = 1; i < genomic_segments.size(); ++i) {
    // adjust [i-1].end [i].start
    size_t prev_end = genomic_segments[i - 1].end;
    size_t cur_start = genomic_segments[i].start;
    size_t mid_point = prev_end + (cur_start - prev_end) / 2;
    genomic_segments[i - 1].end = mid_point;
    genomic_segments[i].start = mid_point;
  }
  return genomic_segments;
}

/**
 * @brief populates the .num_obvs and .mean_logr fields of a list of GenomicSemgents. Also populate num_snps and
 * mean_dh if those fields are present
 * @param genomic_segments
 * @param logrs
 * @return same genomic_segments object, after populating specified fields
 */
std::vector<GenomicSegment>& PopulateGenomicSegmentOptionalFields(std::vector<GenomicSegment>& genomic_segments,
                                                                  const Observations& obvs) {
  // we need to turn logrs into an interval tree
  // look up each genomic segment inside the interval tree
  // then update the mean and num logrs
  IntervalTrees trees(obvs);
  for (auto& gseg : genomic_segments) {
    std::vector<size_t> idxs = trees.LookUp(gseg.contig, gseg.start, gseg.end);
    if (idxs.empty()) {
      Logging::Info("Segment {}:{}-{} not represented in observations! Will not update optional fields...",
                    gseg.contig,
                    gseg.start,
                    gseg.end);
      if (obvs.IsBAFSegObvs()) {
        gseg.num_snps = 0;
        gseg.mean_dh = NAN;
      } else {
        gseg.num_obs = 0;
        gseg.mean_logr = 0;
        gseg.arr_start = *(std::min_element(idxs.begin(), idxs.end()));
        gseg.arr_end = *(std::max_element(idxs.begin(), idxs.end())) + 1;
      }
    } else {
      if (!obvs.IsBAFSegObvs()) {
        // Populate logr related fields
        f64 mean = 0;
        for (size_t i : idxs) {
          mean += obvs.obvs[i];
        }
        mean = mean / static_cast<f64>(idxs.size());
        gseg.num_obs = idxs.size();
        gseg.mean_logr = mean;
      } else {
        // Populate BAF related fields, DP-weighted DH mean
        f64 total_dp = 0;
        f64 af_dp_sum = 0;
        gseg.num_snps = idxs.size();
        for (size_t i : idxs) {
          total_dp += sqrt(obvs.dps[i]);
          af_dp_sum += sqrt(obvs.dps[i]) * obvs.obvs[i];
        }
        if (total_dp > 0) {
          gseg.mean_dh = af_dp_sum / total_dp;
        } else {
          gseg.mean_dh = NAN;
        }
      }
    }
  }
  return genomic_segments;
}

bool IsSortedGenomic(const std::vector<GenomicSegment>& v) {
  std::set<std::string> seen_contigs;
  for (size_t i = 1; i < v.size(); ++i) {
    if (v[i].contig != v[i - 1].contig) {
      if (seen_contigs.find(v[i].contig) != seen_contigs.end()) {
        Logging::Warn("contigs not ordered");
        return false;
      }
      seen_contigs.insert(v[i - 1].contig);
    } else {
      if (v[i].start <= v[i - 1].start) {
        Logging::Warn("starts not ordered");
        return false;
      }
      if (v[i].end <= v[i - 1].end) {
        Logging::Warn("ends not ordered");
        return false;
      }
    }
  }
  return true;
}

std::vector<GenomicSegment>& SortSegments(std::vector<GenomicSegment>& segments,
                                          const std::vector<std::string>& contigs) {
  std::unordered_map<std::string, size_t> contig_to_order;
  for (size_t i = 0; i < contigs.size(); ++i) {
    contig_to_order[contigs[i]] = i;
  }
  std::sort(segments.begin(), segments.end(), [&contig_to_order](const GenomicSegment& lhs, const GenomicSegment& rhs) {
    return std::tie(contig_to_order[lhs.contig], lhs.start, lhs.end) <
           std::tie(contig_to_order[rhs.contig], rhs.start, rhs.end);
  });
  return segments;
}
}  // namespace xoos::cnc::segmentation
