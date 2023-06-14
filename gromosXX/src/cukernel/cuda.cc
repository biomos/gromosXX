/**
 * @file cuda.cc
 * globals of the cukernel library
 */

#include "../stdheader.h"
#include "config.h"

double cukernel_ver = 0.10;

namespace cukernel
{
  char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int kernel_debug_level = 0;
  int constraints_debug_level = 0;
  int interaction_debug_level = 0;
  int pairlist_debug_level = 0;
  int utils_debug_level = 0;
#endif
}