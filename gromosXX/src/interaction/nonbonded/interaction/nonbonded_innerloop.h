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
    /**
     * Constructor
     */
    explicit Nonbonded_Innerloop(Nonbonded_Parameter &nbp) : 
       m_param(&nbp), m_charge_shape(-1) {}
    
    /**
     * accessor to charge shape
     */
    inline void charge_shape(int a) { m_charge_shape = a; }
    inline int charge_shape() { return m_charge_shape; }
    
    /**
     * (normal) interaction
     */
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
     * 1-4 interaction
     * (always shortrange)
     */
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
    void RF_solvent_interaction_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     topology::Chargegroup_Iterator const & cg_it,
     math::Periodicity<t_nonbonded_spec::boundary_type> const &periodicity
     );

    /**
     * fast innerloop for SPC water model
     */
    void spc_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int start,
     int end,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
    /**
     * Calculation of the electric field (polarization)
     */
    void electric_field_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
    );
    
    /**
     * Calculation of the self energy (polarization)
     */
    void self_energy_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
    );
    
    /**
     * calculation of the short range lennard-jones and
     * real space lattice sum electrostatics
     */
    void lj_ls_real_innerloop
    (
     topology::Topology & topo, 
     configuration::Configuration & conf,
     unsigned int i,
     unsigned int j,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
    /**
     * calculation of the real space lattice sum electrostatics
     * for excluded neighbors
     */
    void ls_real_excluded_innerloop
    (
     topology::Topology & topo, 
     configuration::Configuration & conf,
     unsigned int i,
     unsigned int j,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
  protected:
    Nonbonded_Parameter * m_param;
    int m_charge_shape;
  };
} // interaction

#include "nonbonded_innerloop.cc"

#endif
