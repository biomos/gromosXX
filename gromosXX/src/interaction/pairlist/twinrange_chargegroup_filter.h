/**
 * @file twinrange_chargegroup_filter.h
 * filter for twin range using chargegroups (and exclusions).
 */

#ifndef INCLUDED_TWINRANGE_CHARGEGROUP_FILTER_H
#define INCLUDED_TWINRANGE_CHARGEGROUP_FITLER_H

namespace interaction
{
  /**
   * @class Twinrange_Chargegroup_Filter
   * provide filtering using a twinrange cutoff
   * with chargegroups.
   */
  template<typename t_simulation, typename t_base, typename t_innerloop>
  class Twinrange_Chargegroup_Filter
    : public Twinrange_Filter<t_simulation, t_base, t_innerloop>
  {
  public:
    /**
     * Constructor.
     */
    Twinrange_Chargegroup_Filter(t_base &base);

    void prepare(t_simulation &sim);
    
    bool range_chargegroup_pair(t_simulation const &sim, 
				size_t const i,
				size_t const j,
				simulation::chargegroup_iterator const &it_i,
				simulation::chargegroup_iterator const &it_j);
    
  protected:

    math::VArray m_cg_cog;
    
  };
  
} // interaction

#include "twinrange_chargegroup_filter.tcc"

#endif


  
