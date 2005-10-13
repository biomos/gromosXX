/**
 * @file vgrid.h
 * nonbonded grid
 */

#ifndef INCLUDED_VGRID_H
#define INCLUDED_VGRID_H

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}
namespace simulation
{
  class Simulation;
}

namespace interaction
{
  class Nonbonded_Parameter;
  
  void grid_prepare(topology::Topology const & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation const & sim,
		    Nonbonded_Parameter & param);

  void grid_update(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation const & sim);

  extern double grid_timing;

} // interaction

#endif
