#ifndef INCLUDED_BRBR_SHELL_H
#define INCLUDED_BRBR_SHELL_H

#include <cmath>

namespace interaction {
namespace brbr_shell {

struct Weight {
  double value;
  double derivative; // d(value)/dr, where r and r0 are in nm
};

// Same coordination switch as colvar switching_function for m=2*n.
// Evaluate with reciprocal powers outside r0 to avoid overflow. No hard cutoff.
inline Weight weight(double r, double r0, int n) {
  if (r == 0.0) return Weight{1.0, 0.0}; // supported n >= 2
  const double x = r / r0;
  const double p = std::pow(x <= 1.0 ? x : 1.0 / x, n);
  const double den = 1.0 + p;
  return Weight{x <= 1.0 ? 1.0 / den : p / den,
                -n * p / (r * den * den)};
}

struct Correction {
  double delta; // s_ij - 1
  double dvi;   // d(s_ij)/d(v_i)
  double dvj;
};

inline Correction correction(double scale, double vi, double vj) {
  return Correction{(scale - 1.0) * vi * vj,
                    (scale - 1.0) * vj, (scale - 1.0) * vi};
}

} // namespace brbr_shell
} // namespace interaction
#endif
