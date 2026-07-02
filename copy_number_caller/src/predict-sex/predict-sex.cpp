#include "predict-sex/predict-sex.h"

#include <fstream>
#include <string>

#include <xoos/util/file-functions.h>

#include "coverage.h"
#include "sex.h"

namespace xoos::cnc {
s32 PredictSexMain(const PredictSexOptions& options) {
  file::CheckFilePermissionsAndOutputPathExistence({}, {options.predict_sex_out_fname});
  std::ifstream cov_istream(options.coverage_fname);
  CoverageRecords covs(cov_istream);

  XYRatioResult xy_ratio_result = covs.GetXYRatio();
  Sex sex = covs.PredictSex();
  std::string sex_str = SexToStr(sex);

  std::ofstream ofs(options.predict_sex_out_fname);
  ofs << "X_average_coverage" << "\t" << "Y_average_coverage" << "\t" << "XY_ratio" << "\t" << "Sex" << std::endl;
  ofs << xy_ratio_result.chr_x_avg_cov << "\t" << xy_ratio_result.chr_y_avg_cov << "\t" << xy_ratio_result.xy_ratio
      << "\t" << sex_str << std::endl;
  return EXIT_SUCCESS;
}
}  // namespace xoos::cnc
