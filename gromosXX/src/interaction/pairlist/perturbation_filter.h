/**
 * @file perturbation_filter.h
 * filter for perturbation (and construct a pairlist)
 */

#ifndef INCLUDED_PERTURBATION_FILTER_H
#define INCLUDED_PERTURBATION_FILTER_H

namespace interaction
{
  /**
   * @class Perturbation_Filter
   * provide filtering for perturbation
   */
  template<typename t_simulation, typename t_base, 
	   bool do_perturbation = false>
  class Perturbation_Filter
    : public Basic_Filter<t_simulation, t_base>
  {
  public:
    /**
     * Constructor.
     */
    Perturbation_Filter(t_base &base);
    
    /**
     * prepare the filter.
     */
    void prepare(t_simulation &sim);
    /**
     * filter for perturbed atoms.
     */
    bool perturbed_atom(t_simulation const &sim, size_t const i);
    
  protected:
    
  };
  
} // interaction

#include "perturbation_filter.tcc"

#endif

  
