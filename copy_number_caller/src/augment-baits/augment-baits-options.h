#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {

const bool kAugmentBaitsDefaultWholeGenome = false;
const size_t kAugmentBaitsDefaultWholeGenomeIntervalSize = 1000;
const f64 kAugmentBaitsDefaultMinGcContent = 0.15;
const f64 kAugmentBaitsDefaultMaxGcContent = 0.75;
const f64 kAugmentBaitsDefaultWholeGenomeMinMappability = 0.6;
const size_t kAugmentBaitsDefaultOnTargetIntervalSize = 400;
const size_t kAugmentBaitsDefaultOffTargetMinWidth = 20000;
const size_t kAugmentBaitsDefaultOffTargetTrimLength = 500;
const size_t kAugmentBaitsDefaultOffTargetIntervalSize = 200000;
const f64 kAugmentBaitsDefaultOnTargetMinMappability = 0.6;
const f64 kAugmentBaitsDefaultOffTargetMinMappability = 0.1;
const f64 kAugmentBaitsDefaultSexChromosomeTargetMinMappability = 0.9;
const bool kAugmentBaitsDefaultNoOffTargets = false;
const bool kAugmentBaitsDefaultIncludeAltContigs = false;
const size_t kAugmentBaitsGermlineWGSDefaultIntervalSize = 500;
const f64 kAugmentBaitsGermlineWGSDefaultWholeGenomeMinMappabiilty = 0.8;

struct AugmentBaitsOptions {
  size_t whole_genome_interval_size{};
  f64 whole_genome_min_mappability{};
  size_t on_target_interval_size = 400;
  size_t off_target_min_width = 20000;
  size_t off_target_trim_length = 500;
  size_t off_target_interval_size = 200000;
  f64 min_gc_content = kAugmentBaitsDefaultMinGcContent;
  f64 max_gc_content = kAugmentBaitsDefaultMaxGcContent;
  f64 on_target_min_mappability = 0.6;
  f64 off_target_min_mappability = 0.1;
  f64 sex_chromosome_min_mappability = 0.9;
  bool generate_whole_genome_baits = false;
  bool no_off_targets = false;
  bool include_alt_contigs = false;
};
}  // namespace xoos::cnc
