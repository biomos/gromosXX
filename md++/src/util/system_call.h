/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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

