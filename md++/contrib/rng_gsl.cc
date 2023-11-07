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
 * @file rng_gsl.cc
 *
 * this programs helps to find the right flag number for the GSL
 * random number generator.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rng_gsl
 * @section rng_gsl flag number for random number generator.
 * @date 29. 10. 2008
 *
 * This programs helps to find the right flag number for the GSL
 * random number generator.
 * 
 */
#include <iostream>
#include <iomanip>
#include <string>

#include <gsl/gsl_rng.h>


int main() {
  const gsl_rng_type **t, **t0;
          
  t0 = gsl_rng_types_setup ();
          
  std::cout << "Available Random Number Generators in GSL" << std::endl
            << "-----------------------------------------" << std::endl;
         
  bool hasDefault = false; 
  unsigned int i;
  for(i = 0, t = t0; *t != 0; t++, i++) {
    std::string name((*t)->name);
    if (name == "mt19937")
      hasDefault = true;
    std::cout << std::setw(5) << std::left << i 
              << std::setw(20) << std::left << name
              << std::endl;
  }

  if (hasDefault)
    std::cout << std::setw(5) << std::left << -1
              << std::setw(20) << std::left << "mt19937 (default)"
              << std::endl;
  else 
    std::cout << "no default algorithm available!" << std::endl;

  return 0;
}

