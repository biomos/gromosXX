/**
 * @file utils.h
 * utilities
 */

#ifndef INCLUDED_CUKERNEL_UTILS_H
#define INCLUDED_CUKERNEL_UTILS_H
namespace cukernel {
  /**
   * check for an error and print the error message if there was one.
   * @param[in] err_msg the error message
   * @return the error code
   */
  int check_error(const char * err_msg);
}
#endif

