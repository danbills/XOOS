#pragma once

#include <optional>

#include "core/score-calculator.h"
#include "core/variant-id.h"
#include "util/locked-tsv-writer.h"

namespace xoos::svc {

void WriteShapValueTsvHeader(LockedTsvWriter& writer,
                             const vec<std::string>& feature_names,
                             const std::optional<io::CommandLineInfo>& cmd_info);

vec<std::string> AssembleShapValueTsvRow(const VariantId& vid, const PredictionScore& score);

}  // namespace xoos::svc
