#pragma once

#include <xoos/cli/cli.h>
#include <xoos/types/int.h>

#include "copy-number-caller/copy-number-caller-options.h"

namespace xoos::cnc {

// Adds all copy-number-caller subcommands to the CLI app.
void AddSubcommands(cli::AppPtr app, CopyNumberCallerOptionsPtr& options);

void ConfigureSomaticLikelihoodOptions(CopyNumberCallerOptions& options);
void ConfigureGermlineLikelihoodOptions(CopyNumberCallerOptions& options);
void ConfigureSomaticTumorNormalWGSOptions(cli::ConstAppPtr app, CopyNumberCallerOptions& options);
void ConfigureSomaticTumorTargetedEnrichmentOptions(CopyNumberCallerOptions& options);
void ConfigureGermlineNormalWGSOptions(cli::ConstAppPtr app, CopyNumberCallerOptions& options);

/**
 * @brief Applies post-parse fixups: resolves output file paths relative to
 *        output_dir, resets conditional outputs (BAF, BW), and derives sample_id.
 *        Called automatically by CopyNumberCallerCliMain; expose here for tests.
 */
void CopyNumberCallerCliPreCallback(cli::ConstAppPtr app, CopyNumberCallerOptions& options);

s32 CopyNumberCallerCliMain(s32 argc, char** argv);
}  // namespace xoos::cnc
