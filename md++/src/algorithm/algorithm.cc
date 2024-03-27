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
 * @file algorithm.cc
 * globals of the algorithm library
 */

#include "../stdheader.h"
#include "config.h"

double algorithm_ver = 0.10;

namespace algorithm
{
  char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int algorithm_debug_level = 0;
  int constraints_debug_level = 0;
  int integration_debug_level = 0;
  int temperature_debug_level = 0;
  int pressure_debug_level = 0;
  int virtualatoms_debug_level = 0;
#endif
}
