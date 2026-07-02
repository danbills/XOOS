#pragma once

#include <unordered_set>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/util/hash.h>

#include "core/read-collapser-options.h"
#include "core/region.h"
#include "io/alignment-io.h"

namespace xoos::read_collapser {

/// Per-thread set of read names (qnames) of primary alignments that were marked as duplicate
/// and have an SA tag. Used to propagate duplicate status to supplementary alignments.
using DuplicateReadNames = std::unordered_set<std::string, util::hash::TransparentStringHash, std::equal_to<>>;

/**
 * Performed optimized duplicate marking on a BAM file. This mode will:
 *  1. Partition the BAM file into super regions
 *  2. Mark positional duplicates according to start and end position and wiggle room
 *  3. Produce multiple BAM files, one per super region
 *  4. Optionally add cluster ID and cluster size to the BAM files
 *  5. Produce clustering metrics
 *  6. Merge the output BAM files into a single sorted BAM file if requested
 *
 * @param options Configuration for duplicate marking
 * @return The number of super regions processed
 */
u32 MarkDuplicate(const ReadCollapserOptions& options);

/**
 * Mark duplicates on the BAM file specified in @p options, and merge the per-thread output BAM files
 * into a single sorted BAM file if specified in @p options.
 *
 * @param options Configuration for duplicate marking
 */
void MarkDuplicateAndMergeOutput(const ReadCollapserOptions& options);

/**
 * Perform duplicate marking on a super region by reading alignments from @p alignment_reader for each region in
 * @p super_region, clustering them according to configuration defined in @p options, assigning these clusters
 * an id from @p cluster_id_src, and writing the clusters to @p alignment_writer. Use @p super_regions to check if
 * any reads overlap other super regions.
 *
 * When mark_supplementary_alignments is enabled, the qnames of primary alignments that are marked as duplicate
 * and have an SA tag are inserted into @p duplicate_read_names for later supplementary duplicate marking.
 */
void DuplicateMarkSuperRegion(const ReadCollapserOptions& options,
                              u32 super_region_id,
                              const vec<Region>& super_region,
                              const io::AlignmentReader& alignment_reader,
                              const AlignmentWriter& alignment_writer,
                              const vec<vec<Region>>& super_regions,
                              DuplicateReadNames& duplicate_read_names);

}  // namespace xoos::read_collapser
