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
     * (normal) interaction, new version (no virial calculation)
     */
    void lj_crf_innerloop_2
    (
     topology::Topology & topo, 
     unsigned int i,
     unsigned int j,
     const double dist2,
     double &f,
     double &e_lj, double &e_crf
    );

    /**
     * sasa adjustment
     */
    void sasa_calc_innerloop
    (
            topology::Topology & topo,
            configuration::Configuration & conf,
            unsigned int i,
            simulation::Simulation & sim,
            math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
            );

    /**
     * actual adjustment of the sasa
     */
    void calculate_sasa
    (
            topology::Topology & topo,
            configuration::Configuration & conf,
            bool higher,
            const double & pij,
            unsigned int i,
            unsigned int j,
            const topology::sasa_parameter_struct & sasa_param_i,
            math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
            );

    /**
     * sasa force
     */
    void sasa_force_innerloop
    (
            topology::Topology & topo,
            configuration::Configuration & conf,
            unsigned int i,
            simulation::Simulation & sim,
            math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
            );

    /**
     * calculate sasa forces
     */
    void calculate_sasa_forces
    (
            topology::Topology & topo,
            configuration::Configuration & conf,
            bool higher,
            const double pij,
            unsigned int i,
            unsigned int j,
            math::Vec & force_i,
            const double ri_rh2o,
            double dg_i,
            const topology::sasa_parameter_struct & sasa_param_i,
            simulation::Simulation & sim,
            math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
            );

    /**
     * volume term based on sasa expression
     */
    void sasa_volume_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     simulation::Simulation & sim);

    /**
         * switching function for volume calculation
    */
    void sasa_switching_fct
    (
     configuration::Configuration & conf,
     unsigned int i,
     const double amin,
     const double amax,
     const double adiff
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
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );

    /**
     * Lennard-Jones exception
     * (always shortrange)
     */
    void lj_exception_innerloop
    (
            topology::Topology & topo,
            configuration::Configuration & conf,
            topology::lj_exception_struct const & ljex,
            Storage & storage,
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
     Storage & storage,
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
     Storage & storage,
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
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
     unsigned int eps = 0
     );

    /**
     * fast innerloop for SPC water model (new version)
     */
    void spc_innerloop
    (
     double &e_lj,
     double &e_crf,
     double dist6i,
     double f[9],
     double r2[9],
     double r2i[9],
     double ri[9],
     unsigned int eps = 0
     );
    
    /**
     * fast innerloop for SPC water model using tables (shortrange)
     */
    void shortrange_spc_table_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int start,
     int end,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );

    /**
     * fast innerloop for SPC water model using tables (shortrange, new version)
     */    
    void shortrange_spc_table_innerloop
    (
     double &e_lj,
     double &e_crf,
     double f[9],
     double r2[9]
    );
    
    /**
     * fast innerloop for SPC water model using tables (longrange)
     */
    void longrange_spc_table_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     int start,
     int end,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
     );
    
    /**
     * fast innerloop for SPC water model using tables (longrange, new version)
     */
    void longrange_spc_table_innerloop
    (
     double &e_lj,
     double &e_crf,
     double f[9],
     double r2[9]
    );
    
    /**
     * @struct solvent_pair_parameters
     * cached solvent parameters
     */
    struct solvent_pair_parameter {
      double c12, c6, q;
    };
    
    /**
     * innerloop for generic solvents
     */
    void solvent_innerloop
    (
     topology::Topology & topo,
     solvent_pair_parameter * pair_parameter,
     configuration::Configuration & conf,
     const unsigned int num_solvent_atoms,
     const int start,
     const int end,
     Storage & storage,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
     unsigned int eps = 0
     );
    
    /**
     * Calculation of the electric field (polarisation)
     * and repostion the COS till field is self consistent
     */
    void electric_field_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
     math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
    );
    
    /**
     * Calculation of the self energy (polarisation)
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
     * calculation of the short range lennard-jones only
     */
    void lj_innerloop
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
    
    /**
     * accessor to the nonbonded parameters
     */
    Nonbonded_Parameter * param() {
      return m_param;
    }
    
  protected:
    Nonbonded_Parameter * m_param;
    int m_charge_shape;
  
  private:
    /**
     * helper function to calculate the product of charges
     */
    inline static double charge_product(topology::Topology const & topo, 
            unsigned i, unsigned j);
  };
} // interaction

#include "nonbonded_innerloop.cc"

#endif
