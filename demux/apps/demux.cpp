#include "cli/demux.h"

int main(int argc, char** argv) {
  xoos::cli::StandardMainParam<xoos::demux::DemuxAndTrimParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<xoos::demux::DemuxAndTrimParam>(),
      .define_options = xoos::demux::DefineOptions,
      .main = xoos::demux::DemuxAndTrimPipeline,
      .pre_callback = xoos::demux::CreatePreCallback(),
  };
  return xoos::cli::StandardMain(argc, argv, param);
}
