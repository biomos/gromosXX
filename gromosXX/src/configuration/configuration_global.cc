/**
 * @file configuration_global.cc
 * globals of the configuration library
 */

#include "../stdheader.h"

#include "config.h"

namespace configuration
{
  char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int configuration_debug_level;
  int energy_debug_level;
#endif
}
