/**
 * @file nonbonded_innerloop.h
 * inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_NONBONDED_INNERLOOP_H
#define INCLUDED_NONBONDED_INNERLOOP_H

namespace interaction
{
  /**
   * @class Nonbonded_Innerloop
   * standard non bonded inner loop.
   */
  template<typename t_simulation, typename t_nonbonded_spec>
  class Nonbonded_Innerloop
  {
  public:
    /**
     * copy constructor.
     */
    explicit Nonbonded_Innerloop(Nonbonded_Innerloop<t_simulation, t_nonbonded_spec> 
				 const &nil);
    /**
     * Constructor
     */
    explicit Nonbonded_Innerloop(Nonbonded_Base &base);
    
    /**
     * (normal) interaction
     */
    template<typename t_storage>
    void interaction_innerloop(t_simulation const &sim,
			       size_t const i, size_t const j,
			       t_storage &storage);

    /**
     * 1-4 interaction
     * (always shortrange)
     */
    void one_four_interaction_innerloop(t_simulation &sim,
					size_t const i, size_t const j);
    
    /**
     * RF interaction (solute).
     * (always shortrange)
     */
    void RF_excluded_interaction_innerloop(t_simulation &sim,
					   size_t const i);

    /**
     * RF solvent interaction.
     * (always shortrange)
     */
    void RF_solvent_interaction_innerloop(t_simulation &sim,
					  simulation::chargegroup_iterator const & cg_it);
    
 
  protected:
    Nonbonded_Base &m_base;
    
  };
  
} // interaction

#include "nonbonded_innerloop.tcc"

#endif
