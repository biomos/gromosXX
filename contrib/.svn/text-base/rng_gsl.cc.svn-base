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

