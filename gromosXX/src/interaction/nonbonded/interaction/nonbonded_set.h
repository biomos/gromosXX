/**
 * @file nonbonded_set.h
 * the non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * perturbed and non-perturbed.
 */

#ifndef INCLUDED_NONBONDED_SET_H
#define INCLUDED_NONBONDED_SET_H

namespace interaction
{
  
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Nonbonded_Interaction;
  
  /**
   * @class Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Nonbonded_Set : 
    public Nonbonded_Outerloop<t_interaction_spec>,  
    public Perturbed_Nonbonded_Outerloop<t_interaction_spec, 
					 typename t_perturbation_spec::perturbation_details>,
    public Perturbed_Nonbonded_Pair<t_interaction_spec, typename t_perturbation_spec::perturbation_details>
  {
  public:    

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Nonbonded_Set(Nonbonded_Interaction<t_interaction_spec, t_perturbation_spec> & nbi);
    
    /**
     * initialize some things
     */
    void initialize(topology::Topology const & topo,
		    configuration::Configuration const & conf,
		    simulation::Simulation const & sim,
		    bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    int calculate_interactions(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       int tid = 0, int num_threads = 1);

    /**
     * calculate the hessian for a given atom.
     */
    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  size_t const atom_i, size_t const atom_j,
			  math::Matrix & hessian);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     size_t const i, size_t const j,
			     int pc = -1);

    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    size_t const i, size_t const j,
			    Periodicity_type const & periodicity, 
			    int pc = -1);
    Storage & shortrange_storage()
    {
      return m_shortrange_storage;
    }
    Storage & longrange_storage()
    {
      return m_longrange_storage;
    }
    
    Pairlist & pairlist() { return m_pairlist; }
    Pairlist const & pairlist()const { return m_pairlist; }
    Pairlist & perturbed_pairlist() { return m_perturbed_pairlist; }
    Pairlist const & perturbed_pairlist()const { return m_perturbed_pairlist; }
    

  protected:
    Nonbonded_Interaction<
      t_interaction_spec, t_perturbation_spec> * m_nonbonded_interaction;
    
    Storage m_shortrange_storage;
    Storage m_longrange_storage;

    Pairlist m_pairlist;
    Pairlist m_perturbed_pairlist;
  };
  
} // interaction

// template methods
#include "nonbonded_set.cc"

#endif
