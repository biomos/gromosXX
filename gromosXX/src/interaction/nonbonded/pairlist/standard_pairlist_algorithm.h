/**
 * @file standard_pairlist_algorithm.h
 * create an atomic pairlist with a
 * chargegroup based cut-off criterion.
 */

#ifndef INCLUDED_STANDARD_PAIRLIST_ALGORITHM_H
#define INCLUDED_STANDARD_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Standard_Pairlist_Algorithm
   * creates a pairlist.
   */
  template<typename t_nonbonded_spec>
  class Standard_Pairlist_Algorithm : 
    public Pairlist_Algorithm<t_nonbonded_spec>,
    public t_nonbonded_spec::exclusion_filter_type,
    public t_nonbonded_spec::range_filter_type
  {
  public:
    typedef math::Periodicity<t_nonbonded_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Standard_Pairlist_Algorithm();
    
    /**
     * update the pairlist(s).
     */
    template<typename t_nonbonded_interaction>
    void update(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim, 
		t_nonbonded_interaction &nonbonded_interaction);
    
  protected:

    template<typename t_nonbonded_interaction>
    void do_cg_interaction(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   t_nonbonded_interaction &nonbonded_interaction,
			   topology::Chargegroup_Iterator const &cg1,
			   topology::Chargegroup_Iterator const &cg2,
			   Periodicity_type const & periodicity, int const pc = -1);
    
    template<typename t_nonbonded_interaction>
    void do_cg_interaction_excl(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				t_nonbonded_interaction &nonbonded_interaction,
				topology::Chargegroup_Iterator const &cg1,
				topology::Chargegroup_Iterator const &cg2,
				Periodicity_type const & periodicity, int const pc = -1);

    template<typename t_nonbonded_interaction>
    void do_cg_interaction_inv_excl(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    t_nonbonded_interaction &nonbonded_interaction,
				    topology::Chargegroup_Iterator const &cg1,
				    topology::Chargegroup_Iterator const &cg2,
				    Periodicity_type const & periodicity, int const pc = -1);

    template<typename t_nonbonded_interaction>
    void do_cg_interaction_intra(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim,
				 t_nonbonded_interaction &nonbonded_interaction,
				 topology::Chargegroup_Iterator const &cg1,
				 Periodicity_type const & periodicity, int const pc = -1);

  };
} // interaction

#include "standard_pairlist_algorithm.tcc"

#endif
