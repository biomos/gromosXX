/**
 * @file math.cc
 * globals of the math library
 */

#include "config.h"

namespace math
{
  char const id[] = MD_VERSION;

  double h_bar = 0.0635078;

#ifndef NDEBUG
  int debug_level = 0;
#endif
}
