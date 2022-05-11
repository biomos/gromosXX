/**
 * @file peptide_tutorial_eq_test.h
 * Class to perform tests on the equilibration step of the peptide tutorial
 */
#ifndef INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H
#define	INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H

#include "../stdheader.h"
#include "simulation_test.h"

namespace testing {
  
class Peptide_Tutorial_Eq_Test : public Simulation_Test {

public:

  Peptide_Tutorial_Eq_Test(); 

  virtual ~Peptide_Tutorial_Eq_Test() = default;

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

  void check_parameter_init() override;

};


  
}

#endif /* INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H */