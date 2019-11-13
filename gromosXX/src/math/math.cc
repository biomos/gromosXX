/**
 * @file math.cc
 * globals of the math library
 */

#include "../stdheader.h"

#include "config.h"

double math_ver = 0.10;

namespace math
{
  char const id[] = MD_VERSION;
  const char * get_id() { return id; }

  double h_bar = 0.0635078;

  double spd_l = 299792.458;

  double k_Boltzmann = 0.00831441;

  double four_pi_eps_i = 138.9354;
  
  double eps0_i = 1745.9137;

  /**
   * Avogadro constant, Bohr radius and hartree energy
   * IUPAC. Compendium of Chemical Terminology, 2nd ed.
   * (the "Gold Book"). Compiled by A. D. McNaught and A. Wilkinson.
   * Blackwell Scientific Publications, Oxford (1997).
   * Online version (2019-) created by S. J. Chalk. ISBN 0-9678550-9-8.
   * https://doi.org/10.1351/goldbook.
   */
  double avogadro = 6.02214179e23; /* mol-1 */
  double bohr = 5.29177249e-2; /* nm */
  double hartree = 4.3597482e-18; /* J */

#ifndef NDEBUG
  int debug_level = 0;
  int math_debug_level = 0;
#endif
}
