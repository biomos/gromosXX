/**
 * @file system_call.h
 * method for launching an external process
 */
#ifndef SYSTEM_CALL_H
#define	SYSTEM_CALL_H

namespace util {
  extern std::string default_filename_argument;
  /**
   * function that calls an external command, provides standard input from input
   * file and write the output to output file
   * @param command
   * @param input_file the file from which the input is read.
   *        Provide an empty string to skip this.
   * @param output_file the file to which the output is written.
   *        Provide an empty string to let the function create a temporary file and provide its name in the string
   * @return 0 on success, non-zero on failure.
   */
  int system_call(const std::string & command,
          std::string & input_file = default_filename_argument,
          std::string & output_file = default_filename_argument);
}

#endif	/* SYSTEM_CALL_H */

