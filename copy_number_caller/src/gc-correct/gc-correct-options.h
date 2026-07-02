#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {

const f64 kGCCorrectDefaultFirstSpan = 0.03;

struct GCCorrectOptions {
  f64 first_span = kGCCorrectDefaultFirstSpan;
};

}  // namespace xoos::cnc
