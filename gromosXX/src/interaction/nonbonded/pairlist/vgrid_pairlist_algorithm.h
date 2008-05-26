/**
 * @file vgrid_pairlist_algorithm.h
 * vgrid pairlist algorithm (reference implementation)
 */

#ifndef INCLUDED_VGRID_PAIRLIST_ALGORITHM_H
#define INCLUDED_VGRID_PAIRLIST_ALGORITHM_H

namespace interaction
{
  class Pairlist;
   
  /**
   * @class VGrid_Pairlist_Algorithm
   * create a triple - range (atom - range) pairlist
   * (suited for vectorization)
   * using a grid based algorithm
   * and calculate long-range forces
   */
  class VGrid_Pairlist_Algorithm : 
    public Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    VGrid_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~VGrid_Pairlist_Algorithm() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    /**
     * prepare the pairlists
     */    
    virtual int prepare
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation &sim
     );
    
    /**
     * update the pairlist
     */
    virtual void update
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );
    
    /**
     * update the pairlist, separating perturbed and nonperturbed interactions
     */
    virtual void update_perturbed
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     interaction::PairlistContainer & perturbed_pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );
  };
} // interaction

#endif
