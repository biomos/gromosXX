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
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Range_Filter
    : public Filter
  {
  public:
    typedef Chargegroup_Grid<t_interaction_spec::boundary_type> Chargegroup_Grid_type;
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
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
    void prepare_cog(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim);
    /**
     * check the distance of a chargegroup pair.
     * add to the nonbonded_interaction (longrange)
     * if between short- and longrange cutoff,
     * return false if shorter, true if longer.
     */
    bool range_chargegroup_pair
    (topology::Topology & topo, configuration::Configuration & conf,
     simulation::Simulation & sim, 
     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
     size_t const i, size_t const j,
     topology::Chargegroup_Iterator const &it_i,
     topology::Chargegroup_Iterator const &it_j,
     Periodicity_type const & periodicity);

    /**
     * check the distance of a chargegroup pair.
     * add to the nonbonded_interaction (longrange)
     * if between short- and longrange cutoff,
     * return false if shorter, true if longer.
     * shift the first chargegroup by a shift vector
     * instead of calculating the nearest image.
     */
    bool range_chargegroup_pair
    (topology::Topology & topo,	configuration::Configuration & conf,
     simulation::Simulation & sim, 
     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
     size_t const i, size_t const j,
     topology::Chargegroup_Iterator const &it_i,
     topology::Chargegroup_Iterator const &it_j,
     int pc,
     Periodicity_type const & periodicity);

    /**
     * check the distance between two atoms.
     * add to nonbonded_interaction longrange, or filter.
     */
    bool range_atom_pair
    (topology::Topology & topo, configuration::Configuration & conf,
     simulation::Simulation & sim,
     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
     size_t const i, size_t const j,
     Periodicity_type const & periodicity);

    /**
     * check the distance between two atoms.
     * add to nonbonded_interaction longrange, or filter.
     */
    bool range_atom_pair
    (topology::Topology & topo, configuration::Configuration & conf,
     simulation::Simulation & sim,
     Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
     size_t const i, size_t const j, int pc,
     Periodicity_type const & periodicity);
    
    /**
     * put the chargegroup center of geometries
     * on the grid.
     */
    void grid_cog(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  Chargegroup_Grid_type & grid);

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

  private:

  };
  
} // interaction

#include "range_filter.cc"

#endif
