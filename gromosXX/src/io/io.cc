/**
 * @file io.cc
 * globals of the io library
 */

#include "config.h"

namespace io
{
  char const id[] = MD_VERSION;

#ifndef NDEBUG
  int debug_level = 0;
  int trajectory_debug_level = 0;
  int input_debug_level = 0;
  int topology_debug_level = 0;
#endif
}
