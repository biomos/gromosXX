/**
 * @file interaction.cc
 * globals of the interaction library
 */

#include "../stdheader.h"
#include "config.h"

double interaction_ver = 0.10;

namespace interaction
{
  extern char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int forcefield_debug_level = 0;
  int interaction_debug_level = 0;
  int pairlist_debug_level = 0;
  int filter_debug_level = 0;
  int nonbonded_debug_level = 0;
  int latticesum_debug_level = 0;
  int bonded_debug_level = 0;
  int special_debug_level = 0;
#endif
}
