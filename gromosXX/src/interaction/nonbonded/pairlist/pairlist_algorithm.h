/**
 * @file pairlist_algorithm.h
 * base class for the pairlist(s) generating
 * algorithms.
 */

#ifndef INCLUDED_PAIRLIST_ALGORITHM_H
#define INCLUDED_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Pairlist_Algorithm
   * creates a pairlist.
   */
  template<typename t_interaction_spec, bool perturbed>
  class Pairlist_Algorithm:
    public Exclusion_Filter<t_interaction_spec>,
    public Range_Filter<t_interaction_spec, perturbed> 
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm() {}
    /**
     * destructor.
     */
    virtual ~Pairlist_Algorithm() {}
   
    /**
     * prepare the pairlist(s).
     */
    virtual void prepare(){ assert(false); }
    
    /**
     * update the pairlist(s).
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation &sim, 
			Nonbonded_Set<t_interaction_spec, perturbed> &nbs,
			size_t begin, size_t end, size_t stride){ assert(false);}
    
  protected:
  };
} // interaction

#endif
