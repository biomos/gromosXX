/**
 * @file twinrange_filter.h
 * filter for twin range (and exclusions).
 */

#ifndef INCLUDED_TWINRANGE_FILTER_H
#define INCLUDED_TWINRANGE_FITLER_H

namespace interaction
{
  /**
   * @class Twinrange_Filter
   * provide filtering using a twinrange cutoff.
   */
  template<
    typename t_simulation,
    /** interaction parameter */
    typename t_base,
    /** the force-field loop */
    typename t_innerloop,
    /** to enable perturbation */
    typename t_basic_filter = Basic_Filter<t_simulation, t_base> >

  class Twinrange_Filter
    : public t_basic_filter,
      public t_innerloop,
      public Storage
  {
  public:
    /**
     * Constructor.
     */
    Twinrange_Filter(t_base &base);

    /**
     * prepare the filter.
     */
    void prepare(t_simulation &sim);
    
    /**
     * filter a pair for the twin-range.
     * atoms outside the outer range are ignored,
     * if the pair is inside the inner range, false is returned,
     * otherwise the interaction is calculated and the longrange forces
     * and energies stored in the filter.
     */
    bool range_pair(t_simulation const &sim, size_t const i,
		     size_t const j);
    
  protected:
    /**
     * the longrange cutoff.
     */
    double m_cutoff_short_2;
    /**
     * the shortrange cutoff.
     */
    double m_cutoff_long_2;
    
  };
  
} // interaction

#include "twinrange_filter.tcc"

#endif


  
