#pragma once
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/segmentation-options.h"

using std::size_t;

namespace xoos::cnc::segmentation {
void ParentSpecificBinarySegmentationMain(const CopyNumberCallerOptions& options);
std::vector<GenomicSegment> ParentSpecificBinarySegmentation(const Observations& logrs,
                                                             const Observations& dh_vals,
                                                             const std::vector<GenomicSegment>& seed_segments,
                                                             const SegmentationOptions& segmentation_options,
                                                             const SampleMetadataOptions& sample_metadata_options,
                                                             size_t n_threads,
                                                             const fs::path& logr_segments_out,
                                                             const io::CommandLineInfo& command_line_info);
std::vector<GenomicSegment> SegmentByChromosome(const Observations& obvs, const SegmentationOptions& options);
std::vector<GenomicSegment> SegmentBySegments(const std::vector<GenomicSegment>& seed_segments,
                                              const Observations& obvs,
                                              const SegmentationOptions& options,
                                              size_t n_threads);
}  // namespace xoos::cnc::segmentation
