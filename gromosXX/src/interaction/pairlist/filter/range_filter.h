/**
 * @file range_filter.h
 * filter for twin range.
 */

#ifndef INCLUDED_RANGE_FILTER_H
#define INCLUDED_RANGE_FITLER_H

namespace interaction
{
  /**
   * @class Range_Filter
   * provide filtering using a twinrange cutoff.
   * provides methods both for atom and chargegroup
   * based cutoffs.
   */
  template<typename t_simulation, typename t_nonbonded_spec>
  class Range_Filter
    : public Filter<t_simulation, t_nonbonded_spec>
  {
  public:
    /**
     * Constructor.
     */
    Range_Filter();

    void set_cutoff(double const cutoff_short, double const cutoff_long);
    
    void prepare_cog(t_simulation &sim);
    
    template<typename t_nonbonded_interaction>
    bool range_chargegroup_pair(t_simulation & sim, 
				t_nonbonded_interaction &nonbonded_interaction,
				size_t const i,
				size_t const j,
				simulation::chargegroup_iterator const &it_i,
				simulation::chargegroup_iterator const &it_j);
    
    template<typename t_nonbonded_interaction>
    bool range_atom_pair(t_simulation & sim,
			 t_nonbonded_interaction &nonbonded_interaction,
			 size_t const i,
			 size_t const j);
    
  protected:

    /**
     * the chargegroup center of geometries.
     */
    math::VArray m_cg_cog;
    /**
     * squared shortrange cutoff.
     */
    double m_cutoff_short_2;
    /**
     * squared longrange cutoff.
     */
    double m_cutoff_long_2;

  };
  
} // interaction

#include "range_filter.tcc"

#endif
