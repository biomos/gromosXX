/**
 * @file math.cc
 * globals of the math library
 */

#include "config.h"

double math_ver = 0.10;

namespace math
{
  char const id[] = MD_VERSION;

  double h_bar = 0.0635078;

  double k_Boltzmann = 0.00831441;

  double four_pi_eps_i = 138.935;

#ifndef NDEBUG
  int debug_level = 0;
  int math_debug_level = 0;
#endif
}
