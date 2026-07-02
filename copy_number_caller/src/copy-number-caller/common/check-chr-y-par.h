#pragma once

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "baits.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
/**
 * @brief Checks whether the chrY PAR region was masked during alignment by looking for any alignments in the chrY PAR
 * region defined in the seed segments file. If any alignments are found in the chrY PAR region, a warning message is
 * logged since this could significantly affect copy number predictions in allosomes. If the seed segments do not
 * specify any segments in the chrY PAR region, a warning is logged since we cannot check whether the chrY PAR region
 * was masked during alignment.
 * @param bam_file BAM file to check for alignments in the chrY PAR region
 * @param seed_segments Seed segments to check for segments in the chrY PAR region
 * @return true if alignments are found in the chrY PAR region, false otherwise
 */
bool CheckChrYPARUnmasked(const fs::path& bam_file, const vec<segmentation::GenomicSegment>& seed_segments);

/**
 * @brief Removes baits that overlap with the chrY PAR regions defined in the seed segments file if CheckChrYPARUnmasked
 * returns true. If CheckChrYPARUnmasked returns false, the original baits are returned without any modifications.
 * @param bam_file BAM file to check for alignments in the chrY PAR region
 * @param baits Bait records to filter if chrY PAR region is unmasked
 * @param seed_segments Seed segments to check for segments in the chrY PAR region
 * @return filtered bait records with baits that overlap with the chrY PAR regions removed if CheckChrYPARUnmasked
 * returns true, original bait records otherwise
 */
BaitRecords RemovePARIfChrYPARUnmasked(const fs::path& bam_file,
                                       const BaitRecords& baits,
                                       const vec<segmentation::GenomicSegment>& seed_segments);
}  // namespace xoos::cnc
