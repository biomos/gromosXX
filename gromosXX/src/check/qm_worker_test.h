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
#include "test_parameters.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_worker.h"
#include "../interaction/qmmm/qm_zone.h"


namespace testing {

class QM_Worker_Test : public ::testing::Test {

public:

  /**
   * @brief Construct a new qm worker test object
   * 
   * @param parameter 
   * @param results 
   */
  QM_Worker_Test(const Parameter& parameter, const Results& results);

  /**
   * @brief Destroy the qm worker test object
   * 
   */
  virtual ~QM_Worker_Test() = default;

protected:

  /**
   * Called before the execution of each test
   * 
   */
  virtual void SetUp() override = 0;

  /**
   * Called after the execution of each test
   * 
   */
  virtual void TearDown() override = 0;

 /**
   * @brief Initializes the parameters
   * 
   */
  virtual void init_parameters() = 0;

  /**
   * @brief Initializes the expected results
   * 
   */
  virtual void init_results() = 0;

  /**
   * @brief Initializes the expected results for parameters read from the input file
   * 
   */
  virtual void init_results_parameters() = 0;

  /**
   * @brief Initializes the expected results for units and conversion factors
   * 
   */
  virtual void init_results_units() = 0;

  /**
   * @brief Initializes the expected results for element and iac atom mapping
   * 
   */
  virtual void init_results_elements() = 0;

  /**
   * @brief Initializes the expected results for the initial QM zone
   * 
   */
  virtual void init_results_qm_zone_init() = 0;

  /**
   * @brief Checks if Gromos has read in the input parameters correctly
   * 
   */
  virtual void check_parameter_init();

  /**
   * Instance of a test simulation based on the input files
   * 
   */
  Test_Simulation test_sim_;

  /**
   * Structure that stores the expected test results
   * 
   */
  Results results_;

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
