/**
 * @file test_parameters.h
 * define parameters for testing
 */

#ifndef INCLUDED_TEST_PARAMETERS_H
#define INCLUDED_TEST_PARAMETERS_H

#include <cstdint>
#include "../stdheader.h"
#include "check.h"

namespace testing {

/**
 * An enum that stores error codes that may occur during testing
 * 
 */
enum Error {

  /**
   * code for success
   * 
   */
  success,

  /**
   * code for failure
   * 
   */
  failure,

  /**
   * error during parsing of arguments
   * 
   */
  parse_arguments,

  /**
   * error during reading of input files
   * 
   */
  read_input,

  /**
   * error during initialization of simulation
   * 
   */
  md_init,

  /**
   * error during run of simulation
   * 
   */
  md_run,

  /**
   * error caused by too low number of steps provided for simulation
   * 
   */
  md_steps

};

/**
 * A struct that stores parameters for a test simulation
 * 
 */
struct Parameter {

  /**
   * @brief Construct a new Parameter object
   * 
   * @param name Name of the test
   * @param binary_name Name of the binary
   */
  Parameter(const std::string& name, const std::string& binary_name, const std::string& root) 
    : name_(name), binary_name_(binary_name), root_(root) {}

  /**
   * @brief Adds a feature (input, conf, topo, ...) and associated input file to the struct
   * 
   * @param feature Feature (input, conf, topo, ...)
   * @param file Name of the file (relative to root_)
   */
  void add_input(const std::string& feature, const std::string& file) {
    input_.push_back(std::make_pair(feature, TOP_SOURCE_DIR "/" + root_ + "/" + file));
  }

  /**
   * name of the test simulation
   * 
   */
  const std::string name_;

  /**
   * name of the binary
   * 
   */
  std::string binary_name_;

  /**
   * relative path to the root of the test file
   * 
   */
  const std::string root_;

  /**
   * vector of pairs of input features / file names
   * 
   */
  std::vector<std::pair<std::string, std::string>> input_;

};

/**
 * @brief An enums with keys for expected test results
 * 
 */

enum Key {

  kinetic_energy_final,

  potential_energy_final,

  lj_potential_final,

  crf_potential_final,

  qm_energy_init,

  qm_energy_final,

  qmmm,

  qm_ch,

  qm_software,

  qm_cutoff,

  qm_lj,

  qm_constraint,

  qm_mm_scale,

  unit_factor_length,

  unit_factor_energy,

  unit_factor_force,

  unit_factor_charge,

  binary,

  input_file,

  input_coordinate_file,

  input_pointcharges_file,

  output_file,

  output_gradient_file,

  output_mm_gradient_file,

  qm_worker_name,

  qm_zone_size_init,

  mm_zone_size_init,

  qm_zone_charge,

  qm_zone_spin_mult,

  qm_atom_charge_init,

  qm_atom_force_init,

  mm_atom_force_init,

  mm_atom_cos_force_init,

  symbol_h,

  symbol_c,

  symbol_o,

  symbol_p,

  symbol_pd,
  
  element_1,

  element_6,

  element_8,

  element_15,

  element_20,

  element_26,

  element_29,

  element_46

};


/**
 * @brief A struct that stores expected test results
 * 
 */
struct Results {

  /**
   * @brief Results of type std::string
   * 
   */
  std::unordered_map<Key, std::string> strs_;

  /**
   * @brief Results of type double
   * 
   */
  std::unordered_map<Key, double> doubles_;

  /**
   * @brief Results of type int
   * 
   */
  std::unordered_map<Key, int> ints_;

  /**
   * @brief A code for nullptr
   * 
   */
  const nullptr_t UNINITIALIZED = nullptr;

};

/**
 * @brief Stores parameters to test the mechanical embedding scheme with the Orca worker
 * 
 */
extern Parameter orca_parameter_mechanical;

/**
 * @brief Stores expected results for tests of the mechanical embedding scheme with the Orca worker
 * 
 */
extern Results orca_results_mechanical;

/**
 * @brief Stores parameters to test the electrostatic embedding scheme with the Orca worker
 * 
 */
extern Parameter orca_parameter_electrostatic;

/**
 * @brief Stores expected results for tests of the electrostatic embedding scheme with the Orca worker
 * 
 */
extern Results orca_results_electrostatic;

} // namespace testing

#endif /* INCLUDED_TEST_PARAMETERS_H */