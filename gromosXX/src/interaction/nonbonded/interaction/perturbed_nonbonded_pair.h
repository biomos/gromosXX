/**
 * @file perturbed_nonbonded_pair.h
 * perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_PAIR_H
#define INCLUDED_PERTURBED_NONBONDED_PAIR_H

namespace math
{
  template<math::boundary_enum b>
  class Periodicity;
}

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Pair
   * perturbed non bonded pairs.
   */
  class Perturbed_Nonbonded_Pair
  {
  public:
    
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Pair(Nonbonded_Parameter &nbp);

    /**
     * calculate the perturbed pair contributions.
     */
    void perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

  private:
    template<typename t_perturbation_details>
    void _perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_details>
    void _split_perturbed_pair_outerloop(topology::Topology & topo,
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
    template<typename t_interaction_spec, typename t_perturbation_details, math::boundary_enum t_boundary_type>
    void perturbed_pair_interaction_innerloop
    ( topology::Topology & topo, configuration::Configuration & conf,
      simulation::Simulation & sim,
      std::vector<topology::perturbed_two_body_term_struct>
      ::const_iterator const &it,
      math::Periodicity<t_boundary_type> const & periodicity);
 
    Nonbonded_Parameter *m_param;
    Nonbonded_Term m_nonbonded_term;
    Perturbed_Nonbonded_Term m_perturbed_nonbonded_term;

  };
  
} // interaction

#endif
