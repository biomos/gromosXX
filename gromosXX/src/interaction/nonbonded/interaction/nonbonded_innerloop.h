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
  template<typename t_nonbonded_spec>
  class Nonbonded_Innerloop
  {
  public:
    typedef math::Periodicity<t_nonbonded_spec::boundary_type> Periodicity_type;
    
    /**
     * copy constructor.
     */
    explicit Nonbonded_Innerloop(Nonbonded_Innerloop<t_nonbonded_spec> 
				 const &nil);
    /**
     * Constructor
     */
    explicit Nonbonded_Innerloop(Nonbonded_Base &base) : m_base(base) {}
    
    /**
     * (normal) interaction
     */
    template<typename t_storage>
    void interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     size_t const i, size_t const j,
     t_storage & storage,
     Periodicity_type const & periodicity,
     int pc = -1);

    /**
     * 1-4 interaction
     * (always shortrange)
     */
    void one_four_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     size_t const i, size_t const j,
     Periodicity_type const &periodicity);
    
    /**
     * RF interaction (solute).
     * (always shortrange)
     */
    void RF_excluded_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     size_t const i,
     Periodicity_type const & periodicity);

    /**
     * RF solvent interaction.
     * (always shortrange)
     */
    void RF_solvent_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     topology::Chargegroup_Iterator const & cg_it,
     Periodicity_type const &periodicity);
    
 
  protected:
    Nonbonded_Base &m_base;

    /**
     * avoid vector constructor calls?
     */
    // math::Vec r;
    // math::Vec f;
    
  };
  
} // interaction

#include "nonbonded_innerloop.tcc"

#endif
