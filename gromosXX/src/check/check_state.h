/**
 * @file check_state.h
 * a header
 */

namespace check
{
  /**
   * check the state.
   */
  int check_state(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  interaction::Forcefield & ff);
  
} // check

