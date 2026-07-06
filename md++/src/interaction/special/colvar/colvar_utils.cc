/**
 * @file colvar_utils.cc
 * @brief Utility math functions used by collective variables (colvars).
 *
 * This file provides small, performance-oriented helpers used by multiple
 * Colvar implementations (e.g. coordination number, switching functions).
 *
 * Design principles:
 * - Keep functions stateless and reusable
 * - Avoid unnecessary overhead (e.g. std::pow for integer exponents)
 * - Handle numerically delicate regions explicitly (e.g. r ≈ 1)
 *
 * IMPORTANT:
 * All derivatives returned here are with respect to the input argument
 * of the function. If the caller uses a transformed variable (e.g. r/rcut),
 * the chain rule must be applied outside this file.
 */

#include "../../stdheader.h"
#include "colvar_utils.h"
#include <cmath>

namespace interaction {

/**
 * @brief Fast integer power function using exponentiation by squaring.
 *
 * Computes:
 *   base^exp
 *
 * Supports negative exponents by inverting the base.
 *
 * Advantages over std::pow:
 * - Faster for integer exponents
 * - No floating-point exponent overhead
 *
 * @param base Base value
 * @param exp  Integer exponent (can be negative)
 * @return     base^exp
 */
double fast_pow(double base, int exp) {
  if (exp < 0) {
    exp = -exp;
    base = 1.0 / base;
  }

  double result = 1.0;

  while (exp) {
    if (exp & 1)
      result *= base;

    exp >>= 1;
    base *= base;
  }

  return result;
}


/**
 * @brief Smooth switching function commonly used in coordination numbers.
 *
 * Computes:
 *
 *            1 - r^N
 * f(r) = ----------------
 *         1 - r^M
 *
 * where:
 *   r = scaled distance (typically r / r_cut)
 *
 * This function smoothly transitions from ~1 (r → 0) to ~0 (r → ∞),
 * with steepness controlled by parameters N and M.
 *
 * Special case:
 * If M = 2N, the function simplifies to:
 *
 *   f(r) = 1 / (1 + r^N)
 *
 * which is evaluated using a more efficient and numerically stable form.
 *
 * @param rdist  Input value (usually r / r_cut)
 * @param dfunc  Output: derivative df/drdist
 * @param nn     Exponent N
 * @param mm     Exponent M
 *
 * @return       Function value f(rdist)
 *
 * @note
 * The derivative returned is:
 *   dfunc = df / d(rdist)
 *
 * If the caller uses:
 *   rdist = r / r_cut
 *
 * then the physical derivative must be:
 *
 *   df/dr = (dfunc / r_cut)
 *
 * (chain rule must be applied outside this function)
 *
 * @note Numerical stability:
 * Special handling is applied for rdist ≈ 1 to avoid division by zero.
 */
double switching_function(double rdist, double &dfunc, int nn, int mm) {

  // Machine epsilon for double precision
  const double epsilon(std::numeric_limits<double>::epsilon());

  double result;

  /**
   * Special case: M = 2N
   *
   * Simplifies to:
   *   f(r) = 1 / (1 + r^N)
   *
   * This avoids cancellation errors and is faster.
   */
  if (2 * nn == mm) {
    double rNdist = fast_pow(rdist, nn - 1);
    double iden = 1.0 / (1 + rNdist * rdist);

    dfunc = -nn * rNdist * iden * iden;
    result = iden;
  }
  else {

    /**
     * Handle rdist ≈ 1:
     *
     * The expression becomes numerically unstable because:
     *   numerator → 0
     *   denominator → 0
     *
     * Use analytical limit instead.
     */
    if (rdist > (1.0 - 100.0 * epsilon) &&
        rdist < (1.0 + 100.0 * epsilon)) {

      result = static_cast<double>(nn) / static_cast<double>(mm);

      dfunc = 0.5 * nn * (nn - mm) / mm;
    }
    else {

      double rNdist = fast_pow(rdist, nn - 1);
      double rMdist = fast_pow(rdist, mm - 1);

      double num  = 1.0 - rNdist * rdist;   // 1 - r^N
      double iden = 1.0 / (1.0 - rMdist * rdist); // 1 / (1 - r^M)

      double func = num * iden;

      result = func;

      /**
       * Derivative:
       *
       * df/dr = d(num/den) using quotient rule
       *
       * Expanded and simplified form used for efficiency.
       */
      dfunc = ((-nn * rNdist * iden)
             + (func * (iden * mm) * rMdist));
    }
  }

  return result;
}

} // namespace interaction