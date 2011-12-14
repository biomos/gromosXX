/**
 * @file debug.cc
 */

#include "../stdheader.h"

#ifndef NDEBUG

#include "debug.h"

int debug_level = 0;
namespace util
{
  int debug_level = 0;
  int util_debug_level = 0;
  int replica_debug_level = 0;
  int leus_debug_level = 0;
  int bs_leus_debug_level = 0;
}

#endif

double util_ver = 0.10;
