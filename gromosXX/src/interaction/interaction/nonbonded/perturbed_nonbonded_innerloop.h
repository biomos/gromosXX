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
  template<typename t_simulation, typename t_nonbonded_spec>
  class Perturbed_Nonbonded_Innerloop
  {
  public:
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Innerloop(Nonbonded_Base &base);
    
    /**
     * perturbed interaction
     */
    template<typename t_storage>
    void perturbed_interaction_innerloop(t_simulation &sim,
					 size_t const i, size_t const j,
					 t_storage &storage);
    /**
     * perturbed 1-4 interaction
     * (always shortrange)
     */
    void perturbed_one_four_interaction_innerloop(t_simulation &sim,
						  size_t const i, size_t const j);
    
    /**
     * perturbed RF interaction (solute).
     * (always shortrange)
     */
    void perturbed_RF_excluded_interaction_innerloop(t_simulation &sim,
						     std::map<size_t, simulation::Perturbed_Atom>
						     ::const_iterator const & mit);
    /**
     * perturbed pairs.
     * (always shortrange)
     * NO RANGE FILTER FOR PERTURBED PAIRS ??
     * NO SCALING for PERTURBED PAIRS ??
     * NO MOLECULAR VIRIAL CONTRIBUTION ??
     */
    void perturbed_pair_interaction_innerloop(t_simulation &sim,
					      std::vector<simulation::
					      Perturbed_Atompair>::const_iterator const &it);
 
  protected:
    Nonbonded_Base &m_base;
    
  };
  
} // interaction

#include "perturbed_nonbonded_innerloop.tcc"

#endif
