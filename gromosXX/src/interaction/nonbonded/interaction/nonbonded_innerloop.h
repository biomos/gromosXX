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
  class Nonbonded_Innerloop:
    public Nonbonded_Term
  {
  public:
    typedef math::Periodicity<t_nonbonded_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    explicit Nonbonded_Innerloop(Nonbonded_Parameter &nbp) : m_param(&nbp) {}
    
    /**
     * (normal) interaction
     */
    void lj_crf_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Storage & storage,
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
     unsigned int i,
     Periodicity_type const & periodicity);

    /**
     * RF solvent interaction.
     * (always shortrange)
     */
    void RF_solvent_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     topology::Chargegroup_Iterator const & cg_it,
     Periodicity_type const &periodicity);
    
    /**
     * make the hessian available
     */
    using Nonbonded_Term::lj_crf_hessian;
 
  protected:
    Nonbonded_Parameter * m_param;

  };
  
} // interaction

#include "nonbonded_innerloop.cc"

#endif
