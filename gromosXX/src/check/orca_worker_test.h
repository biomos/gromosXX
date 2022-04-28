/**
 * @file orca_worker_test.h
 * Class to perform tests on the Orca QM worker
 */
#ifndef INCLUDED_ORCA_WORKER_T_H
#define	INCLUDED_ORCA_WORKER_T_H

#include "qm_worker_test.h"

namespace testing {

/**
 * @brief A base class for testing the Orca QM worker
 * 
 */
class Orca_Worker_Test : public QM_Worker_Test {

public:

  /**
   * @brief Construct a new Orca_Worker_Test object
   * 
   * @param parameter 
   * @param results 
   */
  Orca_Worker_Test(const Parameter& parameter, const Results& results) : QM_Worker_Test(parameter, results) {}

  /**
   * @brief Destroy the Orca_Worker_Test object
   * 
   */
  virtual ~Orca_Worker_Test() = default;

protected:

  /**
   * @brief Initializes expected results for the file names
   * 
   */
  virtual void init_results_files() = 0;

  /**
   * @brief Initializes expected results for the binary and the worker name
   * 
   */
  virtual void init_results_binary_name();

  /**
   * @brief Initializes expected results for units and conversion factors
   * 
   */
  virtual void init_results_units() override;

  /**
   * @brief Initializes expected results for mapping of element symbols and iac atoms
   * 
   */
  virtual void init_results_elements() override;

  /**
   * @brief Checks if Gromos has read in the input parameters correctly
   * 
   */
  virtual void check_parameter_init() override;

  /**
   * @brief Checks if element symbols and iac atoms are mapped correcly
   * 
   */
  virtual void check_element_mapping();

  /**
   * @brief Checks if the binary and QM worker name have been parsed correcly
   * 
   */
  virtual void check_binary_name();

  /**
   * @brief Checks if input and output files for the QM worker have been parsed correcly
   * 
   */
  virtual void check_files();

  /**
   * @brief Checks if the QM atoms have been initialized correctly
   * 
   */
  virtual void check_qm_atoms_init() = 0;

  /**
   * @brief Checks if the MM atoms have been initialized correcly
   * 
   */
  virtual void check_mm_atoms_init() = 0;

};

} // namespace testing

#endif	/* ORCA_WORKER_T_H */
