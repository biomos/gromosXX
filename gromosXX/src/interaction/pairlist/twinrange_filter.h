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
  template<typename t_simulation, typename t_base, typename t_innerloop>
  class Twinrange_Filter
    : public Basic_Filter<t_simulation>,
      public t_innerloop,
      public Storage
  {
  public:
    /**
     * Constructor.
     */
    Twinrange_Filter(t_base &base);

    void prepare(t_simulation const &sim);
    
    bool range_pair(t_simulation const &sim, size_t const i,
		     size_t const j);
    
  protected:
    t_base &m_base;

    double m_cutoff_short_2;
    double m_cutoff_long_2;
    
  };
  
} // interaction

#include "twinrange_filter.tcc"

#endif


  
