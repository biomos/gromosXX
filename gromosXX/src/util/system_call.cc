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
 * @file system_call.cc
 *
 * implements calls to external programs and handles temporary files creation
 */
#include <cstdlib>
#include <iostream>

#include <unistd.h>

#include "system_call.h"

#include "../../config.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

namespace util
{
  const std::string default_tmp_path = ".";
}

int util::system_call(const std::string & command) {
#ifdef HAVE_SYSTEM
  int system_return = system(command.c_str());
  return WEXITSTATUS(system_return);
#else
  throw std::runtime_error("System call is not implemented for this platform.");
#endif
}

int util::create_tmpfile(std::string & tmp_file) {
#ifdef HAVE_MKSTEMP
  std::string tmp_s;
  if (tmp_file.empty()) {
      tmp_s = util::get_tmppath() + "/gromos-XXXXXX";
  }
  else {
    tmp_s = tmp_file;
  }
  char *tmp_c = &tmp_s[0];
  int fd = mkstemp(tmp_c);
  if (fd >= 1) {
    tmp_file = tmp_s;
  }
  return fd;
#else
  throw std::runtime_error("mkstemp function is not available on this platform.");
#endif
}

int util::create_tmpdir(std::string & tmp_dir) {
#ifdef HAVE_MKDTEMP
  std::string tmp_s;
  if (tmp_dir.empty()) {
    tmp_s = util::get_tmppath() + "/gromos-XXXXXX";
  }
  else {
    tmp_s = tmp_dir;
  }
  char *tmp_c = &tmp_s[0];
  char *dtemp = mkdtemp(tmp_c);
  if(dtemp == nullptr) {
    return 1;
  }
  else {
    tmp_dir = tmp_s;
    return 0;
  }
#else
  throw std::runtime_error("mkdtemp function is not available on this platform.");
#endif
}

std::string util::get_tmppath() {
#ifdef HAVE_GETENV
    char const * tmp_path = getenv("TMPDIR");
    std::string tmp;
    if (tmp_path == nullptr)
      tmp = default_tmp_path;
    else
      tmp = tmp_path;
    return tmp;
#else
  throw std::runtime_error("getenv function is not available on this platform.");
#endif
}
