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

    /**
     * set the short- and long-range cutoff.
     */
    void set_cutoff(double const cutoff_short, double const cutoff_long);
    /**
     * calculate the centers of geometry of the 
     * chargegroups.
     */
    void prepare_cog(t_simulation &sim);
    /**
     * check the distance of a chargegroup pair.
     * add to the nonbonded_interaction (longrange)
     * if between short- and longrange cutoff,
     * return false if shorter, true if longer.
     */
    template<typename t_nonbonded_interaction>
    bool range_chargegroup_pair(t_simulation & sim, 
				t_nonbonded_interaction &nonbonded_interaction,
				size_t const i,
				size_t const j,
				simulation::chargegroup_iterator const &it_i,
				simulation::chargegroup_iterator const &it_j);
    /**
     * check the distance of a chargegroup pair.
     * add to the nonbonded_interaction (longrange)
     * if between short- and longrange cutoff,
     * return false if shorter, true if longer.
     * shift the first chargegroup by a shift vector
     * instead of calculating the nearest image.
     */
    template<typename t_nonbonded_interaction>
    bool range_chargegroup_pair(t_simulation & sim, 
				t_nonbonded_interaction &nonbonded_interaction,
				size_t const i,
				size_t const j,
				simulation::chargegroup_iterator const &it_i,
				simulation::chargegroup_iterator const &it_j,
				int pc);

    /**
     * check the distance between two atoms.
     * add to nonbonded_interaction longrange, or filter.
     */
    template<typename t_nonbonded_interaction>
    bool range_atom_pair(t_simulation & sim,
			 t_nonbonded_interaction &nonbonded_interaction,
			 size_t const i,
			 size_t const j);
    
    void grid_cog(t_simulation const & sim,
		  Chargegroup_Grid<t_simulation> & grid);

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
    /**
     * longrange cutoff.
     */
    double m_cutoff_long;
    /**
     * shortrange cutoff.
     */
    double m_cutoff_short;
    
  };
  
} // interaction

#include "range_filter.tcc"

#endif
