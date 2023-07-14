/**
 * @file utils.h
 * utilities
 */

#ifndef INCLUDED_CUKERNEL_UTILS_H
#define INCLUDED_CUKERNEL_UTILS_H

#include <cstdio>

#define CHECK(call)                                                            \
{   call;                                                                      \
    const cudaError_t error = cudaGetLastError();                              \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}
namespace cukernel {
  /**
   * check for an error and print the error message if there was one.
   * @param[in] err_msg the error message
   * @return the error code
   */
  int check_error(const char * err_msg);
}
#endif

