/**
 * @file simulation.cc
 * globals of the simulation library
 */

#include "config.h"

namespace simulation
{
  char const id[] = MD_VERSION;

#ifndef NDEBUG
  int debug_level;
  int simulation_debug_level;
  int topology_debug_level;
  int system_debug_level;
#endif
}
