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
    // typedef math::Periodicity<t_nonbonded_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    explicit Nonbonded_Innerloop(Nonbonded_Parameter &nbp) : m_param(&nbp) {}
    
    /**
     * (normal) interaction
     */
    // template<typename t_nonbonded_spec>
    void lj_crf_innerloop
    (
     topology::Topology & topo, 
     configuration::Configuration & conf,
     unsigned int i,
     unsigned int j,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );

    /**
     * lennard-jones interaction
     * nearest image free implementation
     * like this only works for molecular virial!
     * this is for interactions in the central computational box
     */
    void lj_crf_innerloop_central
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     unsigned int j,
     Storage & storage
     );

    /**
     * lennard-jones interaction
     * nearest image free implementation
     * like this only works for molecular virial!
     * interactions where one atom has to be shifted
     */
    void lj_crf_innerloop_shift
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     unsigned int j,
     Storage & storage,
     math::Vec const & shift
     );

    /**
     * 1-4 interaction
     * (always shortrange)
     */
    // template<typename t_nonbonded_spec>
    void one_four_interaction_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int i,
     int j,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
    /**
     * RF interaction (solute).
     * (always shortrange)
     */
    // template<typename t_nonbonded_spec>
    void RF_excluded_interaction_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int i,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );

    /**
     * RF solvent interaction.
     * (always shortrange)
     */
    // template<typename t_nonbonded_spec>
    void RF_solvent_interaction_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     topology::Chargegroup_Iterator const & cg_it,
     math::Periodicity<t_nonbonded_spec::boundary_type> const &periodicity
     );

    // template<typename t_nonbonded_spec>
    void spc_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int i,
     int j,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
  protected:
    Nonbonded_Parameter * m_param;
  };
} // interaction

#include "nonbonded_innerloop.cc"

#endif
