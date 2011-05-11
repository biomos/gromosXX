/**
 * @file io.cc
 * globals of the io library
 */

#include "../stdheader.h"
#include "config.h"

double io_ver = 0.10;

namespace io
{
  char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;

  int configuration_debug_level = 0;
  int parameter_debug_level = 0;
  int topology_debug_level = 0;
  int forcefield_debug_level = 0;
#endif
}
