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
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Pairlist_Algorithm:
    public Exclusion_Filter<t_interaction_spec>,
    public Range_Filter<t_interaction_spec, t_perturbation_spec> 
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm() : m_timing(0.0) {}
    /**
     * destructor.
     */
    virtual ~Pairlist_Algorithm() {}
   
    /**
     * prepare the pairlist(s).
     */
    virtual void prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim){ assert(false); }
    
    /**
     * update the pairlist(s).
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation &sim, 
			Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
			unsigned int begin, unsigned int end, unsigned int stride){ assert(false);}
    
    /**
     * timing accessor.
     */
    double timing() { return m_timing; }

  protected:

    double m_timing;

  };
} // interaction

#endif
