/**
 * @file interaction.cc
 * globals of the interaction library
 */

#include "config.h"

namespace interaction
{
  char const id[] = MD_VERSION;

#ifndef NDEBUG
  int debug_level = 0;
  int interaction_debug_level = 0;
#endif
}
