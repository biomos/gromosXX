/**
 * @file prepare_virial.h
 * prepare virial calculation.
 */

#ifndef INCLUDED_PREPARE_VIRIAL_H
#define INCLUDED_PREPARE_VIRIAL_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace util
{
  /**
   * prepare for the virial calculation.
   */
  void prepare_virial(topology::Topology const & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation const & sim);
  
} // util

#endif
