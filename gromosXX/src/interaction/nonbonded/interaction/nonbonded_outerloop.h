/**
 * @file nonbonded_outerloop.h
 * the non bonded outerloops.
 */

#ifndef INCLUDED_NONBONDED_OUTERLOOP_H
#define INCLUDED_NONBONDED_OUTERLOOP_H

namespace interaction
{
  /**
   * @class Nonbonded_Outerloop
   * loops the nonbonded interactions...
   */
  template<typename t_interaction_spec>
  class Nonbonded_Outerloop : 
    public Nonbonded_Innerloop<t_interaction_spec>
  {
  public:    
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Nonbonded_Outerloop(Nonbonded_Parameter & nbp);
    
  protected:

    /**
     * calculate the lj crf interactions.
     */
    void lj_crf_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  std::vector<std::vector<unsigned int> > const & pairlist,
			  Storage & storage);
    /**
     * calculate the 1,4-interactions.
     */
    void one_four_outerloop(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Storage & storage);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    void RF_excluded_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage);
    
  };
  
} // interaction

// template methods
#include "nonbonded_outerloop.cc"

#endif
