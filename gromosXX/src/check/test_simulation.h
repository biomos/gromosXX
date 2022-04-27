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

#include "test_parameters.h"

namespace testing {

class Test_Simulation {

public:
  /**
   * Construct a new test simulation object based on a set of initialization files
   * 
   */
  Test_Simulation(const Parameter& parameter);

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

  /**
   * parameter const accessor
   * 
   * @return const Parameter& 
   */
  const Parameter& parameter() const {
    return parameter_;
  }

  /**
   * parameter accessor
   * 
   * @return Parameter& 
   */
  Parameter& parameter() {
    return parameter_;
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
   * Parameters of the test simulation
   * 
   */
  testing::Parameter parameter_;

};

} // test

#endif	/* TEST_SIMULATION_H */