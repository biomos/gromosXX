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
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Grid_Pairlist_Algorithm : 
    public Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
  {
  public:
    typedef Chargegroup_Grid<t_interaction_spec::boundary_type> Chargegroup_Grid_type;
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor.
     */
    Grid_Pairlist_Algorithm();
    /**
     * Destructor.
     */
    virtual ~Grid_Pairlist_Algorithm(){}
    /**
     * prepare the pairlists
     */
    virtual void prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim);
    
    /**
     * update the pairlist(s).
     */
    void update(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim,
		Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
		unsigned int begin, unsigned int end, unsigned int stride);
    
  private:

    void intra_cell(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
		    std::vector<unsigned int>::const_iterator &cg_st, 
		    std::vector<unsigned int>::const_iterator &cg_to,
		    Periodicity_type const & periodicity);

    template<bool periodic>
    void inter_cell(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
		    std::vector<unsigned int>::const_iterator &cg_st, 
		    std::vector<unsigned int>::const_iterator &cg_to,
		    Chargegroup_Grid_type &grid,
		    int cell[3],
		    Periodicity_type const & periodicity);

  };
} // interaction

#include "grid_pairlist_algorithm.cc"

#endif
