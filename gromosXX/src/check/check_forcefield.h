/**
 * @file check_forcefield.h
 * a header
 */

namespace check
{
  /**
   * check the forcefield.
   */
  int check_forcefield(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       interaction::Forcefield & ff);

  int check_atomic_cutoff(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  interaction::Forcefield & ff);
  
} // check

