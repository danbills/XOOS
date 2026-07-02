#pragma once

#include <xoos/types/fs.h>

namespace xoos::cnc {

struct PredictSexOptions {
  fs::path coverage_fname;
  fs::path predict_sex_out_fname;
};
}  // namespace xoos::cnc
