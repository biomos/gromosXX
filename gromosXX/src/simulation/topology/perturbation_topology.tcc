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
  DEBUG(10, "resizing perturbed_atom " << atoms);
  
  Topology::resize(atoms);
  // standard is non-perturbed atoms
  m_perturbed_atom.resize(atoms, false);
}

inline simulation::Perturbed_Solute &
simulation::Perturbation_Topology::perturbed_solute()
{
  return m_perturbed_solute;
}

inline simulation::Perturbed_Solute const &
simulation::Perturbation_Topology::perturbed_solute()const
{
  return m_perturbed_solute;
}

inline double
simulation::Perturbation_Topology::lambda()
{
  return m_lambda;
}

inline double const 
simulation::Perturbation_Topology::lambda()const
{
  return m_lambda;
}


inline void
simulation::Perturbation_Topology::lambda(double const l)
{
  m_lambda = l;
}
inline double
simulation::Perturbation_Topology::alpha_lj()
{
  return m_alpha_lj;
}

inline void
simulation::Perturbation_Topology::alpha_lj(double const a)
{
  m_alpha_lj = a;
}
inline double
simulation::Perturbation_Topology::alpha_crf()
{
  return m_alpha_crf;
}

inline void
simulation::Perturbation_Topology::alpha_crf(double const a)
{
  m_alpha_crf = a;
}
inline int
simulation::Perturbation_Topology::nlam()
{
  return m_nlam;
}

inline void
simulation::Perturbation_Topology::nlam(int const n)
{
  m_nlam = n;
}

inline std::vector<bool> &
simulation::Perturbation_Topology::perturbed_atom()
{
  return m_perturbed_atom;
}

inline std::vector<bool> const &
simulation::Perturbation_Topology::perturbed_atom()const
{
  return m_perturbed_atom;
}

inline void
simulation::Perturbation_Topology::update_for_lambda()
{
  for(std::map<size_t, simulation::Perturbed_Atom>::const_iterator
	it = perturbed_solute().atoms().begin(),
	to = perturbed_solute().atoms().end();
      it != to; ++it){
    mass()(it->second.sequence_number())
      = (1-lambda()) * it->second.A_mass() + lambda() * it->second.B_mass();
  }
}


