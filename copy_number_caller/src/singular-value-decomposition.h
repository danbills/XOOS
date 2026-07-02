#pragma once

#include <armadillo>
#include <string>

#include <xoos/types/fs.h>

namespace xoos::cnc {
/**
 * @struct construct and hold components for a Singular Value Decomposition of a matrix
 */
struct SingularValueDecomposition {
  SingularValueDecomposition() = default;
  SingularValueDecomposition(const SingularValueDecomposition& rhs) = default;
  SingularValueDecomposition(SingularValueDecomposition&& rhs) noexcept = default;
  SingularValueDecomposition& operator=(const SingularValueDecomposition& rhs) = default;
  SingularValueDecomposition& operator=(SingularValueDecomposition&& rhs) noexcept = default;
  ~SingularValueDecomposition() = default;
  explicit SingularValueDecomposition(const arma::mat& x);
  explicit SingularValueDecomposition(const fs::path& fname, const std::string& group);
  void PerformSingularValueDecomposition(const arma::mat& x);
  void WriteUMatrix(std::ostream& ofs) const;
  void WriteVMatrix(std::ostream& ofs) const;
  void WriteSMatrix(std::ostream& ofs) const;
  void SerializeToHDF5(const fs::path& fname, const std::string& group, bool append) const;
  void LoadFromHDF5(const fs::path& fname, const std::string& group);
  /** @brief truncates u matrix to the first n columsn
   * @param n number of columns [0,n) to keep
   */
  void TruncateU(arma::uword n);
  arma::mat u;
  arma::vec s;
  arma::mat v;
};
}  // namespace xoos::cnc
