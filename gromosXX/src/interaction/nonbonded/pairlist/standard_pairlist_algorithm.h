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
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Standard_Pairlist_Algorithm : 
    public Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
  {
  public:
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Standard_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~Standard_Pairlist_Algorithm(){}
    
    /**
     * prepare the pairlists
     */    
    virtual void prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim);

    /**
     * update the pairlist(s).
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,	
			Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
			size_t begin, size_t end, size_t stride);
        
  protected:

    void do_cg1_loop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
		     topology::Chargegroup_Iterator const & cg1,
		     int cg1_index, int num_solute_cg, int num_cg,
		     Periodicity_type const & periodicity);

    void do_cg_interaction(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
			   topology::Chargegroup_Iterator const &cg1,
			   topology::Chargegroup_Iterator const &cg2,
			   Periodicity_type const & periodicity, int const pc = -1);
    
    void do_cg_interaction_excl(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
				topology::Chargegroup_Iterator const &cg1,
				topology::Chargegroup_Iterator const &cg2,
				Periodicity_type const & periodicity, int const pc = -1);

    void do_cg_interaction_inv_excl(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
				    topology::Chargegroup_Iterator const &cg1,
				    topology::Chargegroup_Iterator const &cg2,
				    Periodicity_type const & periodicity, int const pc = -1);

    void do_cg_interaction_intra(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim,
				 Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
				 topology::Chargegroup_Iterator const &cg1,
				 Periodicity_type const & periodicity, int const pc = -1);


  };
} // interaction

#include "standard_pairlist_algorithm.tcc"

#endif
