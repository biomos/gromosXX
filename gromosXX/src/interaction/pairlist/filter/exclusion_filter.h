/**
 * @file exclusion_filter.h
 * filter for excluded interactions.
 */

#ifndef INCLUDED_EXCLUSION_FILTER_H
#define INCLUDED_EXCLUSION_FILTER_H

namespace interaction
{
  /**
   * @class Exclusion_Filter
   * provide filtering for exclusions.
   */
  template<typename t_simulation, typename t_nonbonded_spec>
  class Exclusion_Filter
    : public Filter<t_simulation, t_nonbonded_spec>
  {
  public:
    /**
     * Constructor.
     */
    Exclusion_Filter();
    /**
     * solute exclusion.
     */
    bool excluded_solute_pair(t_simulation const &sim,
			      size_t const i, size_t const j);

    bool inverse_excluded_solute_pair(t_simulation const &sim,
				      size_t const i, size_t const j);
    
  };
  
} // interaction

#include "exclusion_filter.tcc"

#endif
