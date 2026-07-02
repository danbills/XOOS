#pragma once

#include <optional>
#include <utility>

#include <xoos/types/fs.h>

#include "core/aligner-type.h"
#include "core/config.h"
#include "core/run-type.h"
#include "core/sample-type.h"

namespace xoos::svc {

/**
 * @brief Resolve the tumor-normal-wgs model path from sample type and aligner type.
 * @details Maps the combination of sample type (ffpe, cell-line) and aligner (bwa) to the corresponding
 *          model file in /resources/. Giraffe is not supported for tumor-normal-wgs because no pre-trained
 *          giraffe model exists; pass kCustom with an explicit --model for giraffe-aligned data.
 * @param sample_type The sample preparation type.
 * @param aligner The aligner type. Must not be kGiraffe (unsupported for tumor-normal-wgs) or kCustom
 *                (custom requires explicit --model).
 * @return Path to the model file (e.g., "/resources/model-tumor-normal-wgs-ffpe-bwa.txt.gz").
 * @throws xoos::Error if aligner is kGiraffe or kCustom.
 */
fs::path ResolveTumorNormalModel(SampleType sample_type, AlignerType aligner);

/**
 * @brief Resolve germline SNV and indel model paths from run type and aligner.
 * @details Maps the combination of run type (sbxd, sbxfast) and aligner (bwa, giraffe) to the corresponding
 *          germline model files in /resources/.
 * @param run_type The run type (sbxd or sbxfast).
 * @param aligner The aligner type. Must not be kCustom (custom requires explicit --snv-model/--indel-model).
 * @return Pair of (SNV model path, indel model path).
 * @throws xoos::Error if aligner is kCustom or workflow is not germline.
 */
std::pair<fs::path, fs::path> ResolveGermlineModels(RunType run_type, AlignerType aligner);

/**
 * @brief Resolve germline multi-sample SNV and indel model paths from run type and aligner.
 * @details Maps the combination of run type (sbxd, sbxfast) and aligner (bwa, giraffe) to the corresponding
 *          germline multi-sample model files in /resources/.
 * @param run_type The run type (sbxd or sbxfast).
 * @param aligner The aligner type. Must not be kCustom (custom requires explicit --snv-model/--indel-model).
 * @return Pair of (SNV model path, indel model path).
 * @throws xoos::Error if aligner is kCustom or workflow is not germline.
 */
std::pair<fs::path, fs::path> ResolveGermlineMultiSampleModels(RunType run_type, AlignerType aligner);

// Build a model_thresholds lookup key from sample type and aligner (e.g. "ffpe-bwa").
// Used by ConfigureTumorNormalWgsWorkflow to populate default thresholds and by
// ResolveTumorNormalThresholds to look them up at runtime.
//
// @param sample_type The sample preparation type (e.g. kFfpe, kCellLine).
// @param aligner The aligner type (e.g. kBwa). Must not be kCustom.
// @return A string key in the format "{sample_type}-{aligner}".
std::string TumorNormalModelThresholdsKey(SampleType sample_type, AlignerType aligner);

// Resolve per-model ML score thresholds for tumor-normal-wgs from the config's model_thresholds map.
// Builds a lookup key from sample_type and aligner (e.g. "ffpe-bwa") and searches the map.
// Returns the matching ModelThresholds if found, or std::nullopt if the key is absent or the map is empty.
// @param sample_type The sample preparation type.
// @param aligner The aligner type. Must not be kCustom. Giraffe is rejected upstream by
//                ResolveTumorNormalModel, so only bwa keys are resolved here.
// @param model_thresholds The per-model thresholds map from SVCConfig.
// @return Optional ModelThresholds for the given combination.
std::optional<ModelThresholds> ResolveTumorNormalThresholds(SampleType sample_type,
                                                            AlignerType aligner,
                                                            const StrMap<ModelThresholds>& model_thresholds);

}  // namespace xoos::svc
