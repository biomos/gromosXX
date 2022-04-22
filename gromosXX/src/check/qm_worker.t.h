/**
 * @file qm_worker.t.h
 * Base class to perform tests on QM workers
 */
#ifndef INCLUDED_QM_WORKER_T_H
#define	INCLUDED_QM_WORKER_T_H

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_worker.h"

class QM_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  QM_Worker_Test(std::string binary_name
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file) : binary_name(binary_name)
                                                       , topology_file(topology_file)
                                                       , simulation_file(simulation_file)
                                                       , configuration_file(configuration_file)
                                                       , qmmm_file(qmmm_file)
                                                       , coordinate_trajectory_file(coordinate_trajectory_file)
                                                       , energy_trajectory_file(energy_trajectory_file)
                                                       , final_configuration_file(final_configuration_file) {}

  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~QM_Worker_Test();

  /**
   * Initializes the simulation
   * 
   */
  int init_simulation();

protected:

  /**
   * Instance of the MD simulation
   * 
   */
  algorithm::Algorithm_Sequence md;
  
  /**
   * Instance of the system's topology
   * 
   */
  topology::Topology topo;

  /**
   * Instance of the system's simulation
   * 
   */
  simulation::Simulation sim;

  /**
   * Instance of the system's configuration
   * 
   */
  configuration::Configuration conf;

  /**
   * Store the name of the binary
   * 
   */
  std::string binary_name;

  /**
   * Pointer to the QM/MM interaction
   * 
   */
  interaction::QMMM_Interaction* qmmm_ptr;

  /**
   * Pointer to the QM worker
   * 
   */
  interaction::QM_Worker* qm_worker_ptr;

  /**
   * Path to the topology file
   * 
   */
  std::string topology_file;

  /**
   * Path to the simulation file
   * 
   */
  std::string simulation_file;

  /**
   * Path to the configuration file
   * 
   */
  std::string configuration_file;

  /**
   * Path to the qmmm specification file
   * 
   */
  std::string qmmm_file;

  /**
   * Path to the coordinate trajectory file
   * 
   */
  std::string coordinate_trajectory_file;

  /**
   * Path to the energy trajectory file
   * 
   */
  std::string energy_trajectory_file;

  /**
   * Path to the final configuration file
   * 
   */
  std::string final_configuration_file;

  /**
   * A C++ data structure for command line arguments
   * 
   */
  std::vector<std::string> arguments;

};

#endif	/* QM_WORKER_T_H */
