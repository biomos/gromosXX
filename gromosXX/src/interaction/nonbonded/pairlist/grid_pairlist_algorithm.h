/**
 * @file grid_pairlist_algorithm.h
 * create an atomic pairlist with a
 * chargegroup based cut-off criterion.
 * using a grid!
 */

#ifndef INCLUDED_GRID_PAIRLIST_ALGORITHM_H
#define INCLUDED_GRID_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Grid_Pairlist_Algorithm
   * creates a pairlist.
   */
  template<typename t_nonbonded_spec>
  class Grid_Pairlist_Algorithm : 
    public Standard_Pairlist_Algorithm<t_nonbonded_spec>
  {
  public:
    /**
     * Constructor.
     */
    Grid_Pairlist_Algorithm();
    
    /**
     * update the pairlist(s).
     */
    template<typename t_nonbonded_interaction>
    void update(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim,
		t_nonbonded_interaction &nonbonded_interaction);
    
  private:

    template<typename t_nonbonded_interaction>
    void intra_cell(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    t_nonbonded_interaction & nonbonded_interaction,
		    std::vector<size_t>::const_iterator &cg_st, 
		    std::vector<size_t>::const_iterator &cg_to,
		    math::Periodicity<t_nonbonded_spec::boundary_type> 
		    const & periodicity);

    template<typename t_nonbonded_interaction, bool periodic>
    void inter_cell(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    t_nonbonded_interaction & nonbonded_interaction,
		    std::vector<size_t>::const_iterator &cg_st, 
		    std::vector<size_t>::const_iterator &cg_to,
		    Chargegroup_Grid<t_nonbonded_spec::boundary_type> &grid,
		    int cell[3],
		    math::Periodicity<t_nonbonded_spec::boundary_type> 
		    const & periodicity);

  };
} // interaction

#include "grid_pairlist_algorithm.tcc"

#endif
