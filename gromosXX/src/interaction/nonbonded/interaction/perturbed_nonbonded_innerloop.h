/**
 * @file perturbed_nonbonded_innerloop.h
 * perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_INNERLOOP_H
#define INCLUDED_PERTURBED_NONBONDED_INNERLOOP_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Innerloop
   * perturbed non bonded inner loop.
   */
  template<typename t_interaction_spec,
	   typename t_perturbation_details>
  class Perturbed_Nonbonded_Innerloop:
    public Perturbed_Nonbonded_Term
  {
  public:

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    explicit Perturbed_Nonbonded_Innerloop(Nonbonded_Parameter &nbp): m_param(&nbp){}
    
    /**
     * perturbed interaction
     */
    void perturbed_lj_crf_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Storage &storage,
     Periodicity_type const & periodicity
     );

    /**
     * perturbed 1-4 interaction
     * (always shortrange)
     */
    void perturbed_one_four_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Periodicity_type const & periodicity);
    
    /**
     * perturbed RF interaction (solute).
     * (always shortrange)
     */
    void perturbed_RF_excluded_interaction_innerloop
    ( topology::Topology & topo,
      configuration::Configuration & conf,
      std::map<unsigned int, topology::Perturbed_Atom>::const_iterator const & mit,
      Periodicity_type const & periodicity);
    
    /**
     * Calculation of the perturbed electric field (polarisation)
     */
    void perturbed_electric_field_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
     Periodicity_type const & periodicity
    );
    
    /**
     * Calculation of the perturbed self energy (polarisation)
     */
    void perturbed_self_energy_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i,
     Storage & storage,
     Periodicity_type const & periodicity
    );

  protected:

    Nonbonded_Parameter * m_param;
    
  };
  
} // interaction

#include "perturbed_nonbonded_innerloop.cc"

#endif
