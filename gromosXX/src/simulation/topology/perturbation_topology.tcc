/**
 * @file perturbation_topology.tcc
 * inline methods definition
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE topology

#include "../../debug.h"

/**
 * Constructor
 */
inline simulation::Perturbation_Topology::Perturbation_Topology()
  : Topology(),
    m_perturbed_atom(0)
{
}

inline void simulation::Perturbation_Topology::resize(size_t const atoms)
{
  Topology::resize(atoms);
  // standard is non-perturbed atoms
  m_perturbed_atom.resize(atoms, false);
}

inline simulation::Perturbed_Solute &
simulation::Perturbation_Topology::perturbed_solute()
{
  return m_perturbed_solute;
}

inline double
simulation::Perturbation_Topology::lambda()
{
  return m_lambda;
}

inline void
simulation::Perturbation_Topology::lambda(double const l)
{
  m_lambda = l;
}

