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
