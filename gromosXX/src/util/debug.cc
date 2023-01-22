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
      int timing_debug_level = 0;
      int leus_debug_level = 0;         //Todo: these should be moved out of util! bschroed
      int bs_leus_debug_level = 0;      //Todo: these should be moved out of util! bschroed
    }

#endif

double util_ver = 0.10;
