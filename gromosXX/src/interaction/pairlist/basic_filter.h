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
  template<typename t_simulation>
  class Basic_Filter
  {
  public:
    bool exclusion_solute_pair(t_simulation const &sim, size_t const i,
			       size_t const j);
    bool exclusion_solvent_pair(t_simulation const &sim, size_t const i,
				size_t const j);

  };
  
  
} // interaction

#include "basic_filter.tcc"

#endif
