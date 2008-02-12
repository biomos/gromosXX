/**
 * @file create_simulation.h
 * create a minimum simulation
 */

#ifndef INCLUDED_CREATE_SIMULATION_H
#define INCLUDED_CREATE_SIMULATION_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}
namespace algorithm{
	class Algorithm_Sequence;
}
namespace io{
	class In_Topology;
}

namespace util
{

  /**
   * contains all important entities
   * for an md simulation.
   */
  struct simulation_struct
  {
    topology::Topology topo;
    configuration::Configuration conf;
    simulation::Simulation sim;
    algorithm::Algorithm_Sequence md;
  };
  
  /**
   * create a minimum simulation environment.
   * intended for tests (and maybe Gromos++)
   */
  int create_simulation(std::string topo,
			std::string pttopo,
			std::string conf,
			std::string param,
			util::simulation_struct & sim,
			io::In_Topology & in_topo,
			std::string distrest = "",
			std::string dihrest = "",
			bool quiet = false);
}

#endif
