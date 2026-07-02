#include "core/read-collapser-cli.h"

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/validators/file-extension-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/file-functions.h>

#include "consensus/consensus.h"
#include "core/read-collapser-options.h"
#include "mark-duplicate/mark-duplicate.h"
#include "util/cli-option-util.h"

namespace xoos::read_collapser {

using cli::AddEnumOption;
using cli::AddOptionalEnumOption;
using enum_util::FormatEnumName;
using enum_util::FormatEnumNames;
using util::container::Contains;

void SetCommandLineInfo(const cli::ConstAppPtr app, const ReadCollapserOptionsPtr& options) {
  options->command_line_info = cli::GetCommandLineInfo(app);
}

void DefineConsensusOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options) {
  AddPresetOption(app, kConsensusPresetsMap, kOptGroupNamePresetOptions);
  AddCommonInputOptions(app, options, kOptGroupNameInputOptions);
  AddCommonOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddConsensusOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddCommonReadFilteringOptions(app, options, kOptGroupNameReadFilterOptions);
  AddConsensusReadFilterOptions(app, options, kOptGroupNameReadFilterOptions);
  AddCommonClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddConsensusClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddConsensusOptions(app, options, kOptGroupNameConsensusOptions);
  AddConsensusDebugOptions(app, options, kOptGroupNameConsensusDebugOptions);
  AddCommonPerformanceOptions(app, options, kOptGroupNamePerformanceOptions);
  AddCommonHiddenOptions(app, options);
}

void DefineMarkdupOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options) {
  AddPresetOption(app, kMarkdupPresetsMap, kOptGroupNamePresetOptions);
  AddCommonInputOptions(app, options, kOptGroupNameInputOptions);
  AddCommonOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddMarkdupOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddCommonReadFilteringOptions(app, options, kOptGroupNameReadFilterOptions);
  AddMarkdupReadFilterOptions(app, options, kOptGroupNameReadFilterOptions);
  AddCommonClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddCommonPerformanceOptions(app, options, kOptGroupNamePerformanceOptions);
  AddCommonHiddenOptions(app, options);
}

void DefineOptions(cli::AppPtr app, ReadCollapserOptionsPtr& options) {
  // Add 'markdup' subcommand
  cli::AddSubcommand<ReadCollapserOptions>(app,
                                           "markdup",
                                           DefineMarkdupOptions,
                                           options,
                                           MarkDuplicateAndMergeOutput,
                                           "Identifies, marks, and optionally removes duplicate reads.",
                                           [](const cli::ConstAppPtr app, const ReadCollapserOptionsPtr& options) {
                                             ValidateMarkdupOptions(app, options);
                                             SetCommandLineInfo(app, options);
                                           });

  // Add 'consensus' subcommand
  cli::AddSubcommand<ReadCollapserOptions>(app,
                                           "consensus",
                                           DefineConsensusOptions,
                                           options,
                                           FastClusterAndConsensus,
                                           "Generates high-quality consensus sequences.",
                                           [](const cli::ConstAppPtr app, const ReadCollapserOptionsPtr& options) {
                                             ValidateConsensusOptions(app, options);
                                             SetCommandLineInfo(app, options);
                                           });

  // Require exactly one subcommand to be executed
  app->require_subcommand(1);
}

}  // namespace xoos::read_collapser
