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
    m_perturbed_atom(0),
    m_lambda(0),
    m_nlam(1)
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

inline int const
simulation::Perturbation_Topology::nlam()const
{
  return m_nlam;
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


/**
 * calculate constraint degrees of freedom.
 */
inline void
simulation::Perturbation_Topology::
calculate_constraint_dof(simulation::Multibath &multibath)const
{
  // do nonperturbed constraints
  simulation::Topology::calculate_constraint_dof(multibath);
  
  DEBUG(7, "and the perturbd distance constraints (DOF calc)");

  // substract perturbed constraints
  std::vector<Perturbed_Solute::perturbed_distance_constraint_struct>
    ::const_iterator 
    c_it = perturbed_solute().distance_constraints().begin(),
    c_to = perturbed_solute().distance_constraints().end();
  
  size_t com_bath_i, ir_bath_i, com_bath_j, ir_bath_j;

  for( ; c_it != c_to; ++c_it){
    
    DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
    multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
    multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);

    multibath[ir_bath_i].dof -= 0.5;
    multibath[ir_bath_j].dof -= 0.5;

    multibath[ir_bath_i].ir_dof -= 0.5;
    multibath[ir_bath_j].ir_dof -= 0.5;

    multibath[ir_bath_i].solute_constr_dof += 0.5;
    multibath[ir_bath_j].solute_constr_dof += 0.5;

  }
  
  for(size_t i=0; i<multibath.size(); ++i){
    DEBUG(7, "dof           " << multibath[i].dof);
    DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
  }
  
}
/**
 * interaction matrix for scaled interactions with energy groups
 */
inline std::map<std::pair<int, int>, std::pair<double, double> > const &
simulation::Perturbation_Topology::energy_group_scaling()const
{
  return m_energy_group_scaling;
}
/**
 * interaction matrix for scaled interactions with energy groups
 */
inline std::map<std::pair<int, int>, std::pair<double, double> >  &
simulation::Perturbation_Topology::energy_group_scaling()
{
  return m_energy_group_scaling;
}

/**
 * add solvent molecules to the simulation (system).
 */
inline void simulation::Perturbation_Topology
::solvate(size_t solv, size_t num_molecules)
{

  Topology::solvate(solv, num_molecules);
  m_perturbed_atom.resize(num_atoms(), false);
}
