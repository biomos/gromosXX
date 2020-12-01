/**
 * @file system_call.h
 * method for launching an external process
 */
#ifndef SYSTEM_CALL_H
#define	SYSTEM_CALL_H

namespace util {
  extern const std::string default_tmp_path;
  /**
   * function that calls an external command
   * @param command complete command including optional stdout and stderr pipe
   * @return 0 on success, non-zero on failure.
   */
  int system_call(const std::string & command);

  /**
   * function to generate temporary file
   * requires filename reference as input, overwrites it with the generated filename value
   * and returns file descriptor. Last six characters of string should be XXXXXX. If empty
   * string is given, the function gets path from util::get_tmppath() and uses default name
   * gromos-XXXXXX.
   */
  int create_tmpfile(std::string & tmp_file);
  
  /**
   * function to generate temporary directory
   * requires directory name reference as input and overwrites it with the generated dirname
   * value and returns 0 on success. Last six characters of string should be XXXXXX. If empty
   * string is given, the function gets path from util::get_tmppath() and uses default name
   * gromos-XXXXXX.
   */
  int create_tmpdir(std::string & tmp_dir);

    /**
   * function to get system temporary directory path. If TMPDIR is defined, it returns TMPDIR,
   * otherwise it returns default_tmp_path. TMPDIR should be set to shared memory filesystem
   * (e.g. /dev/shm) to avoid slow I/O.
   */
  std::string get_tmppath();
}

#endif	/* SYSTEM_CALL_H */

