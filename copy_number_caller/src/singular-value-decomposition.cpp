#include "singular-value-decomposition.h"

#include <stdexcept>

#include <xoos/log/logging.h>

namespace xoos::cnc {

SingularValueDecomposition::SingularValueDecomposition(const arma::mat& x) {
  PerformSingularValueDecomposition(x);
}

SingularValueDecomposition::SingularValueDecomposition(const fs::path& fname, const std::string& group) {
  LoadFromHDF5(fname, group);
}

void SingularValueDecomposition::PerformSingularValueDecomposition(const arma::mat& x) {
  if (!arma::svd(u, s, v, x)) {
    throw(std::runtime_error("SVD failed\n"));
  }
}

void SingularValueDecomposition::TruncateU(arma::uword n) {
  if (n < u.n_cols) {
    u = u.cols(0, n - 1);
  } else {
    Logging::Warn("cannot truncate U to {} columns because U only has {} columns", n, u.n_cols);
  }
}

void SingularValueDecomposition::WriteUMatrix(std::ostream& ofs) const {
  u.each_row([&ofs](const arma::rowvec& b) {
    size_t i = 0;
    b.for_each([&ofs, &i](const arma::mat::elem_type& x) {
      if (i) {
        ofs << " " << x;
      } else {
        ofs << x;
      }
      i += 1;
    });
    ofs << std::endl;
  });
}

void SingularValueDecomposition::WriteVMatrix(std::ostream& ofs) const {
  v.each_row([&ofs](const arma::rowvec& b) {
    size_t i = 0;
    b.for_each([&ofs, &i](const arma::mat::elem_type& x) {
      if (i) {
        ofs << " " << x;
      } else {
        ofs << x;
      }
      i += 1;
    });
    ofs << std::endl;
  });
}

void SingularValueDecomposition::WriteSMatrix(std::ostream& ofs) const {
  s.for_each([&ofs](const arma::mat::elem_type& x) { ofs << x << std::endl; });
}

void SingularValueDecomposition::SerializeToHDF5(const fs::path& fname, const std::string& group, bool append) const {
  if (append) {
    u.save(arma::hdf5_name(fname.string(), group + "/svd_u", arma::hdf5_opts::append));
    v.save(arma::hdf5_name(fname.string(), group + "/svd_v", arma::hdf5_opts::append));
    s.save(arma::hdf5_name(fname.string(), group + "/svd_s", arma::hdf5_opts::append));
  } else {
    u.save(arma::hdf5_name(fname.string(), group + "/svd_u", arma::hdf5_opts::replace));
    v.save(arma::hdf5_name(fname.string(), group + "/svd_v", arma::hdf5_opts::append));
    s.save(arma::hdf5_name(fname.string(), group + "/svd_s", arma::hdf5_opts::append));
  }
}

void SingularValueDecomposition::LoadFromHDF5(const fs::path& fname, const std::string& group) {
  u.load(arma::hdf5_name(fname.string(), group + "/svd_u"));
  v.load(arma::hdf5_name(fname.string(), group + "/svd_v"));
  s.load(arma::hdf5_name(fname.string(), group + "/svd_s"));
}

}  // namespace xoos::cnc
