#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <xoos/types/float.h>

#include "observations.h"
#include "sex.h"

namespace xoos::cnc::segmentation {
struct Segment {
  size_t start = 0;  // 0-based inclusive
  size_t end = 0;    // 0-based exclusive
  f64 t = 0;         // t-value for this segment compared to parent
  f64 p = 0;         // p-value for this segment
  f64 n = 0;         // number of values within this segment
  f64 mean = 0;      // mean value in this segment.
};

struct GenomicSegment {
  GenomicSegment() = default;

  GenomicSegment(std::string c, size_t s, size_t e) : contig(std::move(c)), start(s), end(e) {
  }

  GenomicSegment(const Segment& segment, const Observations& obvs);

  void PopulateOptionalFields(const Observations& obvs);
  void SetField(const std::string& field, const std::string& val);
  void SetFieldOrThrowError(const std::string& field, const std::string& val);

  std::vector<std::string> GetPopulatedOptionalCols() const;
  std::string GetFieldAsString(const std::string& field) const;
  std::string contig{};
  size_t start{};  // 0-based inclusive
  size_t end{};    // 0-based exclusive
  std::optional<std::string> id{};
  std::optional<size_t> arr_start{};
  std::optional<size_t> arr_end{};
  std::optional<size_t> num_obs{};
  std::optional<f64> mean_logr{};
  std::optional<size_t> num_snps{};
  std::optional<f64> mean_dh{};
  std::optional<f64> mbaf{};
  std::optional<f64> avg_mean_mapq{};
  std::optional<f64> ploidy{};
  std::optional<f64> purity{};
  std::optional<size_t> total_copy_number{};
  std::optional<f64> expected_total_copy_number{};
  std::optional<size_t> minor_copy_number{};
  std::optional<size_t> major_copy_number{};
  std::optional<f64> logr_likelihood{};
  std::optional<f64> baf_likelihood{};
  std::optional<f64> joint_likelihood{};
  std::optional<Sex> sex{};
  std::optional<std::vector<std::string>> flags{};
  std::optional<bool> in_allosome{};
  std::optional<bool> in_pseudo_autosomal_region{};
};

std::vector<GenomicSegment> ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                             const std::vector<Segment>& segments);
std::vector<GenomicSegment>& ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                              const std::vector<Segment>& segments,
                                                              std::vector<GenomicSegment>& genomic_segment);
std::vector<GenomicSegment>& PopulateGenomicSegmentOptionalFields(std::vector<GenomicSegment>& genomic_segments,
                                                                  const Observations& obvs);
std::vector<GenomicSegment>& AdjustGenomicSegmentBreakpoints(std::vector<GenomicSegment>& genomic_segments);
bool IsSortedGenomic(const std::vector<GenomicSegment>& v);
std::vector<GenomicSegment>& SortSegments(std::vector<GenomicSegment>& segments,
                                          const std::vector<std::string>& contigs);
}  // namespace xoos::cnc::segmentation
