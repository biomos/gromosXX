/**
 * @file perturbation_filter.h
 * filter for perturbed atoms.
 */

#ifndef INCLUDED_PERTURBATION_FILTER_H
#define INCLUDED_PERTURBATION_FILTER_H

namespace interaction
{
  /**
   * @class Perturbation_Filter
   * provide filtering for perturbed atoms.
   */
  template<typename t_nonbonded_spec>
  class Perturbation_Filter
    : public Filter<t_nonbonded_spec>
  {
  public:
    /**
     * Constructor.
     */
    Perturbation_Filter();
    /**
     * solute exclusion.
     */
    bool perturbed_atom(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim, size_t const i);
  };
  
} // interaction

#include "perturbation_filter.cc"

#endif
