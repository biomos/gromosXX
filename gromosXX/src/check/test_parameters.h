/**
 * @file test_parameters.h
 * define parameters for testing
 */

#ifndef INCLUDED_TEST_PARAMETERS_H
#define INCLUDED_TEST_PARAMETERS_H

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

} // namespace testing

#endif /* INCLUDED_TEST_PARAMETERS_H */