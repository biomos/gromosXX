/**
 * @file configuration.cc
 * globals of the configuration library
 */

#include "config.h"

namespace configuration
{
  char const id[] = MD_VERSION;

#ifndef NDEBUG
  int debug_level = 0;
  int configuration_debug_level;
  int energy_debug_level;
#endif
}
