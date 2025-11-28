/**
 * @file colvar_utils.h
 * @brief Utility math functions shared by collective variables (colvars)
 */

#ifndef INCLUDED_COLVAR_UTILS_H
#define INCLUDED_COLVAR_UTILS_H

namespace interaction {

  /**
   * @brief Smooth switching function used for coordination number colvars.
   * 
   * Computes a continuous function f(r) and its derivative df/dr, depending on 
   * the distance ratio r/rcut and exponents n and m.
   *
   * f(r) = (1 - (r/rcut)^n) / (1 - (r/rcut)^m)
   *
   * @param rdist Ratio of distance to cutoff (r/rcut)
   * @param dfunc Output derivative of f with respect to r
   * @param nn Exponent n
   * @param mm Exponent m
   * @return double Value of switching function
   */
  double switching_function(double rdist, double &dfunc, int nn, int mm);

  /**
   * @brief Fast integer exponentiation
   *
   * @param base Base value
   * @param exp Integer exponent
   * @return double Result of base^exp
   */
  double fast_pow(double base, int exp);

} // namespace interaction

#endif
