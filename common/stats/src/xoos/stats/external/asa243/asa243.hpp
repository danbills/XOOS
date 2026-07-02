#pragma once

/**
 * Code originally taken from https://people.sc.fsu.edu/~jburkardt/cpp_src/asa243/asa243.html
 * Original author: J Bukardt, who adapted it from the FORTRAN code presented in the paper:
 * Lenth R. Algorithm AS 243: Cumulative Distribution Function of the Non-Central T Distribution
 * This code was modified by Taher Mun on May 9 2024 and March 12 2026. Changes include
 * 2. putting the code under the xoos::stats namespace (May 2024)
 * 3. adding `NOLINT` checks to keep the function names (March 2026)
 * 4. removing the ifault pointer argument, replacing ifault with runtime errors (May 2024)
 * 5. various changes to support our linters (May 2024)
 */

namespace xoos::stats {
double alnorm(double x, bool upper);                       // NOLINT(readability-identifier-naming)
double betain(double x, double p, double q, double beta);  // NOLINT(readability-identifier-naming)
double tnc(double t, double df, double delta);             // NOLINT(readability-identifier-naming)
}  // namespace xoos::stats
