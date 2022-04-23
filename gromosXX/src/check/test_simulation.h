/**
 * @file test_simulation.h
 * A test simulation that mimics md.cc and can be used for testing
 */
#ifndef INCLUDED_TEST_SIMULATION_H
#define	INCLUDED_TEST_SIMULATION_H

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../io/configuration/out_configuration.h"

namespace test {

class Test_Simulation {

public:
  /**
   * Construct a new test simulation object based on a set of initialization files
   * 
   */
  Test_Simulation(std::string binary_name
               , std::string test_title
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file) : binary_name_(binary_name)
                                                       , traj_(test_title)
                                                       , topology_file_(topology_file)
                                                       , simulation_file_(simulation_file)
                                                       , configuration_file_(configuration_file)
                                                       , qmmm_file_(qmmm_file)
                                                       , coordinate_trajectory_file_(coordinate_trajectory_file)
                                                       , energy_trajectory_file_(energy_trajectory_file)
                                                       , final_configuration_file_(final_configuration_file) {}

                                                       

  /**
   * Destroy the test simulation object
   * 
   */
  virtual ~Test_Simulation() = default;

  /**
   * Initializes the simulation in analogy to md.cc
   * 
   * @return int error code
   */
  int init_simulation();

  /**
   * Runs a single step of the simulation
   * 
   * @return int error code
   */
  int run_single_step();

  /**
   * Runs all steps of the simulation
   * 
   * @return int error code
   */
  int run_simulation();

  /**
   * algorithm_sequence const accessor
   * 
   * @return const algorithm::Algorithm_Sequence& 
   */
  const algorithm::Algorithm_Sequence& md() const {
    return md_;
  }

  /**
   * algorithm_sequence accessor
   * 
   * @return algorithm::Algorithm_Sequence& 
   */
  algorithm::Algorithm_Sequence& md() {
    return md_;
  }
  
  /**
   * topology const accessor
   * 
   * @return const topology::Topology& 
   */
  const topology::Topology& topo() const {
    return topo_;
  }

  /**
   * topology accessor
   * 
   * @return topology::Topology& 
   */
  topology::Topology& topo() {
    return topo_;
  }

  /**
   * simulation const accessor
   * 
   * @return const simulation::Simulation& 
   */
  const simulation::Simulation& sim() const {
    return sim_;
  }

  /**
   * simulation accessor
   * 
   * @return simulation::Simulation& 
   */
  simulation::Simulation& sim() {
    return sim_;
  }

  /**
   * configuration const accessor
   * 
   * @return const configuration::Configuration& 
   */
  const configuration::Configuration& conf() const {
    return conf_;
  }

  /**
   * configuration accessor
   * 
   * @return configuration::Configuration& 
   */
  configuration::Configuration& conf() {
    return conf_;
  }

protected:

  /**
   * Instance of the MD simulation
   * 
   */
  algorithm::Algorithm_Sequence md_;
  
  /**
   * Instance of the system's topology
   * 
   */
  topology::Topology topo_;

  /**
   * Instance of the system's configuration
   * 
   */
  configuration::Configuration conf_;

  /**
   * Instance of the system's simulation
   * 
   */
  simulation::Simulation sim_;

  /**
   * Instance of the system's trajectory
   * 
   */
  io::Out_Configuration traj_;

  /**
   * Store the name of the binary
   * 
   */
  std::string binary_name_;

  /**
   * Path to the topology file
   * 
   */
  std::string topology_file_;

  /**
   * Path to the simulation file
   * 
   */
  std::string simulation_file_;

  /**
   * Path to the configuration file
   * 
   */
  std::string configuration_file_;

  /**
   * Path to the qmmm specification file
   * 
   */
  std::string qmmm_file_;

  /**
   * Path to the coordinate trajectory file
   * 
   */
  std::string coordinate_trajectory_file_;

  /**
   * Path to the energy trajectory file
   * 
   */
  std::string energy_trajectory_file_;

  /**
   * Path to the final configuration file
   * 
   */
  std::string final_configuration_file_;

  /**
   * A C++ data structure for command line arguments
   * 
   */
  std::vector<std::string> arguments_;

};

} // test

#endif	/* TEST_SIMULATION_H */