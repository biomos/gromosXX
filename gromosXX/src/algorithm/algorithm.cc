/**
 * @file algorithm.cc
 * globals of the algorithm library
 */

#include "config.h"

namespace algorithm
{
  char const id[] = MD_VERSION;

#ifndef NDEBUG
  int debug_level = 0;
  int constraint_debug_level = 0;
  int integration_debug_level = 0;
  int algorithm_debug_level = 0;
  int temperature_debug_level = 0;
  
#endif
}
