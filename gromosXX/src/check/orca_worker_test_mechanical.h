/**
 * @file orca_worker_test_mechanical.h
 * Class to perform tests on the Orca QM worker with mechanical embedding
 */
#ifndef INCLUDED_ORCA_WORKER_TEST_MECHCANICAL_H
#define	INCLUDED_ORCA_WORKER_TEST_MECHCANICAL_H

#include "orca_worker_test.h"

namespace testing {

/**
 * @brief A class that tests the Orca_Worker using an mechanical embedding
 * 
 */
class Orca_Worker_Test_Mechanical : public Orca_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Test_Mechanical();


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Test_Mechanical() = default;

protected:

  /**
   * @brief Runs before each test
   * 
   */
  void SetUp() override;

  /**
   * @brief Runs after each test
   * 
   */
  void TearDown() override;

  /**
   * @brief Initializes the parameters for test test simulation
   * 
   */
  void init_parameters() override;

  /**
   * @brief Initializes the expected results for a test simulation
   * 
   */
  void init_results() override;

  /**
   * @brief Initializes the expected results for parameters read from the input file
   * 
   */
  void init_results_parameters() override;

  /**
   * @brief Initializes expected results for the file names
   * 
   */
  void init_results_files() override;

  /**
   * @brief Initializes the expected results for the initial QM zone
   * 
   */
  void init_results_qm_zone_init() override;

};

} // namespace testing

#endif /* INCLUDED_ORCA_WORKER_TEST_MECHCANICAL_H */
  