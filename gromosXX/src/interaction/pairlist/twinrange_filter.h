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

    void prepare(t_simulation &sim);
    
    bool range_pair(t_simulation const &sim, size_t const i,
		     size_t const j);
    
  protected:
    double m_cutoff_short_2;
    double m_cutoff_long_2;
    
  };
  
} // interaction

#include "twinrange_filter.tcc"

#endif


  
