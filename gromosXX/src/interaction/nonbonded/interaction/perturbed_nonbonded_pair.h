/**
 * @file perturbed_nonbonded_pair.h
 * perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_PAIR_H
#define INCLUDED_PERTURBED_NONBONDED_PAIR_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Pair
   * perturbed non bonded pairs.
   */
  template<typename t_interaction_spec, typename perturbation_details>
  class Perturbed_Nonbonded_Pair
  {
  public:
    
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Pair(Nonbonded_Parameter &nbp,
			     Nonbonded_Term & nbt,
			     Perturbed_Nonbonded_Term &pnbt);
    

    /**
     * calculate the perturbed pair contributions.
     */
    void perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

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
    Nonbonded_Parameter *m_param;
    Nonbonded_Term * m_nonbonded_term;
    Perturbed_Nonbonded_Term * m_perturbed_nonbonded_term;

  };
  
} // interaction

#include "perturbed_nonbonded_pair.tcc"

#endif
