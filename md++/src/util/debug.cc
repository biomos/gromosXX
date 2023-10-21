/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
