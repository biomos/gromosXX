/**
 * @file perturbed_shake.tcc
 * contains methods for Perturbed_Shake
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraint

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation>
algorithm::Perturbed_Shake<t_simulation>
::Perturbed_Shake(double const tolerance, int const max_iterations)
  : Shake<t_simulation>::Shake(tolerance, max_iterations)
{
}

/**
 * Destructor.
 */
template<typename t_simulation>
algorithm::Perturbed_Shake<t_simulation>
::~Perturbed_Shake()
{
}

/**
 * shake solute
 */
template<typename t_simulation>
int algorithm::Perturbed_Shake<t_simulation>
::solute(typename simulation_type::topology_type &topo,
	 typename simulation_type::system_type &sys,
	 double dt)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...
  
  // check whether we shake
  if (!topo.solute().distance_constraints().size()) return 0;
  
  sys.constraint_force() = 0.0;

  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  int num_iterations = 0;

  int first = 0;

  skip_now.assign(topo.solute().num_atoms(), false);
  skip_next.assign(topo.solute().num_atoms(), true);
  
  bool convergence = false;
  bool pert_convergence = false;
  
  while(!(convergence && pert_convergence)){

    // non-perturbed constraints
    DEBUG(7, "non-perturbed-constraints");
    convergence = _shake(topo, sys, first, skip_now, skip_next,
			 topo.solute().distance_constraints(), true);

    // and perturbed constraints
    DEBUG(7, "perturbed-constraints");
    pert_convergence = _shake(topo, sys, first, skip_now, skip_next,
			      topo.perturbed_solute().distance_constraints(), 
			      true, 
			      topo.solute().distance_constraints().size());

    if(++num_iterations > max_iterations){
      io::messages.add("SHAKE error. too many iterations",
		       "Shake::solute",
		       io::message::critical);
      throw std::runtime_error("SHAKE failure in solute");
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?

  //constraint_force and free energy derivatives
  size_t k = 0;
  for( ; k < topo.solute().distance_constraints().size(); ++k){

    sys.constraint_force()(k) /=  (dt * dt);
    DEBUG(5, "constraint_force " << 
	  sqrt(dot(sys.constraint_force()(k),
		   sys.constraint_force()(k)) ));

  }

  std::vector<simulation::Perturbed_Solute::
    perturbed_distance_constraint_struct>::const_iterator
    it = topo.perturbed_solute().distance_constraints().begin(),
    to = topo.perturbed_solute().distance_constraints().end();
  
  for(; it != to; ++it, ++k){
    
    sys.constraint_force()(k) /= (dt * dt);
    DEBUG(5, "perturbed constraint_force " << 
	  sqrt(dot(sys.constraint_force()(k),
		   sys.constraint_force()(k)) ));

    sys.lambda_energies().constraint_energy -= 
      sqrt(dot(sys.constraint_force()(k),
	       sys.constraint_force()(k))) *
      (it->B_b0 - it->A_b0);
  }
  
  // shaken velocity
  sys.vel() = (sys.pos() - sys.old_pos()) / dt;

  DEBUG(8, "perturbed shake iterations " << num_iterations);

  return num_iterations;

} // solute

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
template<typename t_simulation>
inline void 
algorithm::Perturbed_Shake<t_simulation>
::add_bond_length_constraints(typename t_simulation::topology_type &topo)
{
  DEBUG(7, "perturbed_shake: add_bond_length_constraints");
  
  Shake<t_simulation>::add_bond_length_constraints(topo);
  
  simulation::Perturbed_Solute & solute = topo.perturbed_solute();
  
  // std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    DEBUG(8, "adding pert constraint " << it->i << " - " << it->j);
    solute.add_distance_constraint(it->i, it->j, m_bond_parameter[it->type].r0,
				   m_bond_parameter[it->B_type].r0);
  }
  solute.bonds().clear();
}

/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Perturbed_Shake<t_simulation>
::add_bond_length_constraints(int iac,
			      std::vector<int> const &atom_iac,
			      typename t_simulation::topology_type &topo)
{
  DEBUG(7, "perturbed_shake: add_bond_length_constraints (iac)");
  Shake<t_simulation>::add_bond_length_constraints(iac, atom_iac, topo);
  
  simulation::Perturbed_Solute & solute = topo.perturbed_solute();
  
  std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_iac[it->i] == iac || atom_iac[it->j] == iac){
      DEBUG(8, "adding pert constraint " << it->i << " - " << it->j);
      solute.add_distance_constraint(it->i, it->j, 
				     m_bond_parameter[it->type].r0,
				     m_bond_parameter[it->B_type].r0);
    }
    else
      bonds.push_back(simulation::Perturbed_Bond(*it));
  }
  solute.bonds() = bonds;
}
    
/**
 * add bonds connecting an atom of mass mass to the
 * constraint vector and remove from the bond vector.
 */
template<typename t_simulation>
inline void
algorithm::Perturbed_Shake<t_simulation>
::add_bond_length_constraints(double mass,
			      math::SArray const &atom_mass,
			      typename t_simulation::topology_type &topo)
{
  DEBUG(7, "perturbed_shake: add_bond_length_constraints (mass)");
  Shake<t_simulation>::add_bond_length_constraints(mass, atom_mass, topo);
  
  simulation::Perturbed_Solute &solute = topo.perturbed_solute();

  std::vector<simulation::Perturbed_Bond> bonds;
  std::vector<simulation::Perturbed_Bond>::iterator
    it = solute.bonds().begin(),
    to = solute.bonds().end();

  for( ; it != to; ++it){
    if(atom_mass(it->i) == mass || atom_mass(it->j) == mass){
      DEBUG(8, "adding pert constraint " << it->i << " - " << it->j);      
      solute.add_distance_constraint(it->i, it->j,
				     m_bond_parameter[it->type].r0,
				     m_bond_parameter[it->B_type].r0);
    }
    else
      bonds.push_back(simulation::Perturbed_Bond(*it));
  }
  solute.bonds() = bonds;
}

template<typename t_simulation>
inline void
algorithm::Perturbed_Shake<t_simulation>
::init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
       io::InInput &input)
{
  Shake<t_simulation>::init(sim, args, topo, input);

  sim.topology().perturbed_solute().
    set_distance_constraints(sim.topology().lambda());

  // give the constraint force the correct size...
  sim.system().constraint_force().resize(sim.topology().solute().
					 distance_constraints().size() + 
					 sim.topology().perturbed_solute().
					 distance_constraints().size());

}
