/**
 * @file debug.cc
 */

#ifndef NDEBUG

#include <util/stdheader.h>

#include "debug.h"

int debug_level = 0;
namespace util
{
  int debug_level = 0;
  int util_debug_level = 0;
}

#endif
