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
  class Nonbonded_Outerloop
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Outerloop(Nonbonded_Parameter & nbp);
    
    /**
     * calculate the lj crf interactions.
     */
    void lj_crf_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Pairlist const & pairlist,
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

    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  unsigned int atom_i, unsigned int atom_j,
			  math::Matrix & hessian,
			  Pairlist const & pairlist);

  private:
    /**
     * the nonbonded parameter.
     */
    Nonbonded_Parameter & m_param;

    template<typename t_interaction_spec>
    void _lj_crf_outerloop(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Pairlist const & pairlist,
			   Storage & storage);

    template<typename t_interaction_spec>
    void _one_four_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Storage & storage);

    template<typename t_interaction_spec>
    void _RF_excluded_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage);

    template<typename t_interaction_spec>
    int _calculate_hessian(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   unsigned int atom_i, unsigned int atom_j,
			   math::Matrix & hessian,
			   Pairlist const & pairlist);
    
  };
  
} // interaction

#endif
