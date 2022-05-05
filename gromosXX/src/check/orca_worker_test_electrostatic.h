/**
 * @file orca_worker_test_electrostatic.h
 * Class to perform tests on the Orca QM worker with electrostatic embedding
 */
#ifndef INCLUDED_ORCA_WORKER_TEST_ELECTROSTATIC_H
#define	INCLUDED_ORCA_WORKER_TEST_ELECTROSTATIC_H

#include "orca_worker_test.h"

namespace testing {

/**
 * @brief A class that tests the Orca_Worker using an electrostatic embedding
 * 
 */
class Orca_Worker_Test_Electrostatic : public Orca_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Test_Electrostatic();


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Test_Electrostatic() = default;

protected:

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
   * @brief Initializes the expected results for energies calculated
   * 
   */
  void init_results_energies() override;

  /**
   * @brief Initializes the expected results for forces calculated
   * 
   */
  void init_results_forces() override;

  /**
   * @brief Initializes the expected results for velocities calculated
   * 
   */
  void init_results_velocities() override;

  /**
   * @brief Initializes the expected results for positions calculated
   * 
   */
  void init_results_positions() override;

  void init_results_baths() override;

  void init_results_bonded_terms() override;

  void init_results_nonbonded_terms() override;

  void init_results_special_terms() override;

  void init_results_mass() override;

  void init_results_temperature() override;

  void init_results_volume() override;

  void init_results_pressure() override;

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

  /**
   * @brief Checks if the QM atoms have been initialized correctly
   * 
   */
  void check_qm_atoms_init() override;

  /**
   * @brief Checks if the MM atoms have been initialized correctly
   * 
   */
  void check_mm_atoms_init() override;

};

}

#endif /* INCLUDED_ORCA_WORKER_TEST_ELECTROSTATIC_H */