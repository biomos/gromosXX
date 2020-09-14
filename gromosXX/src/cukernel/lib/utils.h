/**
 * @file utils.h
 * utilities
 */

#ifndef CUKERNEL_UTILS
#define CUKERNEL_UTILS
namespace cudakernel {
  /**
   * check for an error and print the error message if there was one.
   * @param[in] err_msg the error message
   * @return the error code
   */
  int checkError(const char * err_msg);
}
#endif

