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

  /**
   * recover molecular virial from atomic virial
   */
  void atomic_to_molecular_virial(topology::Topology const & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation const & sim);
  
  /**
   * calculate centre of mass, centre of mass kinetic energy
   * of the (sub) molecules
   */
  void centre_of_mass(topology::Topology const & topo,
		      configuration::Configuration & conf,
		      std::vector<math::Vec> & com_pos,
		      std::vector<math::Matrix> & com_ekin,
          	      simulation::Simulation const & sim);

} // util

#endif
