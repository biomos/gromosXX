/**
 * @file perturbed_solute.tcc
 * inline methods of Perturbed_Solute
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE topology

#include "../../debug.h"

inline std::vector<simulation::Perturbed_Bond> &
simulation::Perturbed_Solute::bonds()
{
  return m_bond;
}
