/**
 * @file basic_filter.h
 * filter the interactions while
 * generating the pairlist.
 */

#ifndef INCLUDED_BASIC_FILTER_H
#define INCLUDED_BASIC_FILTER_H

namespace interaction
{
  /**
   * @class Basic_Filter
   * provide filtering for exclusions.
   */
  template<typename t_simulation, typename t_base>
  class Basic_Filter
  {
  public:
    /**
     * Constructor.
     */
    Basic_Filter(t_base &base);

    /**
     * prepare the filter.
     */
    void prepare(t_simulation &sim);
    
    /**
     * solute exclusion.
     */
    bool exclusion_solute_pair(t_simulation const &sim, size_t const i,
			       size_t const j);
    /**
     * solvent exclusion.
     */
    bool exclusion_solvent_pair(t_simulation const &sim, size_t const i,
				size_t const j);

    /**
     * no perturbation here.
     */
    bool perturbed_atom(t_simulation const &sim, size_t const i) {return false; }

  protected:
    t_base &m_base;
    
  };
  
  
} // interaction

#include "basic_filter.tcc"

#endif
