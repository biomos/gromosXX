/**
 * @file error_codes.h
 * define error values for testing
 */

#ifndef INCLUDED_ERROR_CODES_H
#define INCLUDED_ERROR_CODES_H

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

}

#endif /* ERROR_CODES_H */