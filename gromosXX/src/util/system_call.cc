/**
 * @file system_call.cc
 *
 * implements calls to external programs and handles temporary files creation
 */
#include <cstdlib>
#include <iostream>

#include "system_call.h"

#include "../../config.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

int util::system_call(const std::string & command,
          std::string & input_file,
          std::string & output_file) {
  std::string command_to_launch = command;
  if (!input_file.empty()) {
    command_to_launch += " < " + input_file;
  }
  command_to_launch += " 1> " + output_file + " 2>&1 ";
#ifdef HAVE_SYSTEM
  int system_return = system(command_to_launch.c_str());
  return WEXITSTATUS(system_return);
#else
  throw std::runtime_error("System call is not implemented for this platform.");
#endif
}

int util::create_tmpfile(std::string & tmp_file) {
#ifdef HAVE_MKSTEMP
  std::string tmp_s;
  if (tmp_file.empty()) {
    std::string tmp_path = getenv("TMPDIR");
    if (tmp_path.empty()) {
      tmp_s = default_tmp_path + "/gromos-XXXXXX";
    }
    else {
      tmp_s = tmp_path + "/gromos-XXXXXX";
    }
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
    std::string tmp_path = getenv("TMPDIR");
    if (tmp_path.empty()) {
      tmp_s = default_tmp_path + "/gromos-XXXXXX";
    }
    else {
      tmp_s = tmp_path + "/gromos-XXXXXX";
    }
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
