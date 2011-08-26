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
  int read_input
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   std::ostream & os = std::cout,
   bool quiet = false
   );
  
  /**
   * read the rest of the input files necessary for the repex simulation.
   * calls
   * - read_topology
   * - read_special
   * - read_configuration
   */
  int read_input_repex
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in simulation parameter
   */
  int read_parameter
  (
   io::Argument const &args,
   simulation::Simulation & sim,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in topology, perturbation topology, MD sequence
   */
  int read_topology
  (
   io::Argument const &args,
   topology::Topology &topo,
   simulation::Simulation & sim,
   algorithm::Algorithm_Sequence &md_seq,
   std::ostream & os = std::cout,
   bool quiet = false
   );

  /**
   * read in a configuration
   */
  int read_configuration
  (
   io::Argument const &args,
   topology::Topology &topo,
   configuration::Configuration &conf,
   simulation::Simulation & sim,
   std::ostream & os = std::cout,
   bool quiet = false
   );
}

#endif
