/**
 * @file virtual_grain.h
 * update positions / forces of virtual atoms
 */

#ifndef INCLUDED_VIRTUAL_GRAIN_H
#define INCLUDED_VIRTUAL_GRAIN_H

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}

namespace util
{
  /**
   * update positions of virtual atoms based upon
   * real atom positions
   */
  void update_virtual_pos(topology::Topology & cg_topo,
			  configuration::Configuration & cg_conf,
			  topology::Topology & topo,
			  configuration::Configuration & conf);
  /**
   * distribute the forces from virtual atoms on
   * the real atoms
   */
  void update_virtual_force(topology::Topology & cg_topo,
			    configuration::Configuration & cg_conf,
			    topology::Topology & topo,
			    configuration::Configuration & conf);
}

#endif
