/**
 * @file algorithm.h
 * gathers common include directives for algorithm
 */

#include "integration/leap_frog.h"
#include "constraint/shake.h"

#ifndef NDEBUG
namespace algorithm
{
  extern int debug_level;
  extern int constraint_debug_level;
  extern int integration_debug_level;
}
#endif
