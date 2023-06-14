/**
 * @file utils.cu
 * implementation of utils
 */

#include <cuda.h>
#include <iostream>
#include "utils.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE utils

int cukernel::check_error(const char * err_msg) {
#ifndef NDEBUG
  cudaDeviceSynchronize();
#endif
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    std::cout << "CUDA-ERROR " << err_msg << ": " << cudaGetErrorString(error) << std::endl;
  return (int) error;
}

