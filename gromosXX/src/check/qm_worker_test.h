/**
 * @file qm_worker_test.h
 * Base class to perform tests on QM workers
 */
#ifndef INCLUDED_QM_WORKER_T_H
#define	INCLUDED_QM_WORKER_T_H

#include <gtest/gtest.h>
#include <iostream>

#include "../stdheader.h"

#include "test_simulation.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_worker.h"
#include "../interaction/qmmm/qm_zone.h"


namespace testing {

class QM_Worker_Test : public ::testing::Test {

public:

  QM_Worker_Test(std::string binary_name
               , std::string test_title
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file);

  virtual ~QM_Worker_Test() = default;

protected:

  /**
   * Called before the execution of each test
   * 
   */
  void SetUp() override;

  /**
   * Called after the execution of each test
   * 
   */
  void TearDown() override;

  /**
   * Instance of a test simulation based on the input files
   * 
   */
  Test_Simulation test_sim_;

  /**
   * Pointer to the qmmm_interaction object of the test simulation
   * 
   */
  interaction::QMMM_Interaction* qmmm_interaction_ptr;

  /**
   * Pointer to the qm_worker object of the test simulation
   * 
   */
  interaction::QM_Worker* qm_worker_ptr;

  /**
   * Pointer to the qm_zone object of the test simulation
   * 
   */
  interaction::QM_Zone* qm_zone_ptr;

};

}

#endif	/* QM_WORKER_T_H */
