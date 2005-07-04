/**
 * @file read_input.h
 * input parameters
 */

#ifndef INCLUDED_READ_INPUT_H
#define INCLUDED_READ_INPUT_H

namespace util
{
  struct Replica_Data;
}

namespace io
{
  /**
   * read all input files necessary for the simulation.
   * calls
   * - read_parameter
   * - read_topology
   * - read_special
   * - read_configuration
   */
  int read_input(io::Argument const &args,
		 topology::Topology &topo,
		 configuration::Configuration &conf,
		 simulation::Simulation & sim,
		 algorithm::Algorithm_Sequence &md_seq,
		 std::ostream & os = std::cout);

  /**
   * read replica information
   */
  int read_replica_input
  (
   io::Argument const & args,
   topology::Topology & topo,
   std::vector<configuration::Configuration> & conf,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence & md_seq,
   std::vector<util::Replica_Data> & replica_data,
   std::ostream & os);

  /**
   * read in simulation parameter
   */
  int read_parameter(io::Argument const &args,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout);

  /**
   * read in topology, perturbation topology, MD sequence
   */
  int read_topology(io::Argument const &args,
		    topology::Topology &topo,
		    simulation::Simulation & sim,
		    algorithm::Algorithm_Sequence &md_seq,
		    std::ostream & os = std::cout);

  /**
   * read in a configuration
   */
  int read_configuration(io::Argument const &args,
			 topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation & sim,
			 std::ostream & os = std::cout);

  /**
   * read in a configuration
   */
  int read_replica_configuration
  (
   io::Argument const &args,
   topology::Topology &topo,
   std::vector<configuration::Configuration> & conf,
   simulation::Simulation & sim,
   std::vector<util::Replica_Data> & replica_data,
   std::ostream & os = std::cout
   );
  
}

#endif
