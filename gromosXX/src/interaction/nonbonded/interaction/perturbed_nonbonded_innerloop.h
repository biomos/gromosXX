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
  template<typename t_interaction_spec>
  class Perturbed_Nonbonded_Innerloop
  {
  public:

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Innerloop(Nonbonded_Base &base);
    
    /**
     * perturbed interaction
     */
    template<typename t_storage>
    void perturbed_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     size_t const i, size_t const j, t_storage &storage,
     Periodicity_type const & periodicity,
     int pc = -1);

    /**
     * perturbed 1-4 interaction
     * (always shortrange)
     */
    void perturbed_one_four_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     size_t const i, size_t const j,
     Periodicity_type const & periodicity);
    
    /**
     * perturbed RF interaction (solute).
     * (always shortrange)
     */
    void perturbed_RF_excluded_interaction_innerloop
    ( topology::Topology & topo,
      configuration::Configuration & conf,
      std::map<size_t, topology::Perturbed_Atom>::const_iterator const & mit,
      Periodicity_type const & periodicity);
    /**
     * perturbed pairs.
     * (always shortrange)
     * NO RANGE FILTER FOR PERTURBED PAIRS ??
     * NO SCALING for PERTURBED PAIRS ??
     * NO MOLECULAR VIRIAL CONTRIBUTION ??
     */
    void perturbed_pair_interaction_innerloop
    ( topology::Topology & topo, configuration::Configuration & conf,
      simulation::Simulation & sim,
      std::vector<topology::perturbed_two_body_term_struct>
      ::const_iterator const &it,
      Periodicity_type const & periodicity);
 
  protected:
    Nonbonded_Base &m_base;
    
  };
  
} // interaction

#include "perturbed_nonbonded_innerloop.tcc"

#endif
