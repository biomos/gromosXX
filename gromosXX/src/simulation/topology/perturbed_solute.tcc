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

inline std::vector<simulation::Perturbed_Bond> const &
simulation::Perturbed_Solute::bonds()const
{
  return m_bond;
}

inline std::vector<simulation::Perturbed_Angle> &
simulation::Perturbed_Solute::angles()
{
  return m_angle;
}

inline std::vector<simulation::Perturbed_Angle> const &
simulation::Perturbed_Solute::angles()const
{
  return m_angle;
}

inline std::map<size_t, simulation::Perturbed_Atom> &
simulation::Perturbed_Solute::atoms()
{
  return m_atom;
}

inline std::map<size_t, simulation::Perturbed_Atom> const &
simulation::Perturbed_Solute::atoms()const
{
  return m_atom;
}
