/**
 * @file simple_filter.h
 * filter the interactions while
 * generating the pairlists.
 */

#ifndef SIMPLE_FILTER_H
#define SIMPLE_FILTER_H

namespace interaction
{
  /**
   * @class Simple_Filter
   * provide filtering for exclusions.
   * no use of the concept of chargegroups.
   */
  template<typename t_simulation>
  class Simple_Filter
  {
  public:
    /**
     * filter pairs of solute atoms.
     */
    bool solute_pair(t_simulation const & sim, size_t const i,
		     size_t const j);

    bool solvent_pair(t_simulation const & sim, size_t const i,
		      size_t const j);

  protected:
    /**
     * filters a pair for exclusions (and 1,4 interactions).
     * this only works for solvent if the exclusions are added
     * to the arrays (which will normally not be the case).
     * so: do not call the filter for solvent.
     */
    bool filter_excluded(t_simulation const & sim, size_t const i,
			 size_t const j);
    /**
     * filters a pair for solvent exclusions.
     * checks whether i and j belong to the same solvent molecule.
     */
    bool filter_solvent_excluded(t_simulation const & sim, size_t const i,
				 size_t const j);
    
  };
  
} // interaction

#include "simple_filter.tcc"

#endif
