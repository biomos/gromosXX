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
  template<typename t_simulation, typename t_nonbonded_spec>
  class Standard_Pairlist_Algorithm : 
    public Pairlist_Algorithm<t_simulation, t_nonbonded_spec>,
    public t_nonbonded_spec::exclusion_filter_type,
    public t_nonbonded_spec::range_filter_type
  {
  public:
    /**
     * Constructor.
     */
    Standard_Pairlist_Algorithm();
    
    /**
     * update the pairlist(s).
     */
    template<typename t_nonbonded_interaction>
    void update(t_simulation &sim, t_nonbonded_interaction &nonbonded_interaction);
    
  protected:

    template<typename t_nonbonded_interaction>
    void do_cg_interaction(t_simulation & sim,
			   t_nonbonded_interaction &nonbonded_interaction,
			   simulation::chargegroup_iterator const &cg1,
			   simulation::chargegroup_iterator const &cg2);
    
    template<typename t_nonbonded_interaction>
    void do_cg_interaction_excl(t_simulation & sim,
				t_nonbonded_interaction &nonbonded_interaction,
				simulation::chargegroup_iterator const &cg1,
				simulation::chargegroup_iterator const &cg2);

    template<typename t_nonbonded_interaction>
    void do_cg_interaction_inv_excl(t_simulation & sim,
				    t_nonbonded_interaction &nonbonded_interaction,
				    simulation::chargegroup_iterator const &cg1,
				    simulation::chargegroup_iterator const &cg2);

    template<typename t_nonbonded_interaction>
    void do_cg_interaction_intra(t_simulation & sim,
				 t_nonbonded_interaction &nonbonded_interaction,
				 simulation::chargegroup_iterator const &cg1);

  };
} // interaction

#include "standard_pairlist_algorithm.tcc"

#endif
