/**
 * @file utils.cu
 * implementation of utils
 */

#include <cuda.h>
#include <iostream>
#include "utils.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE utils

int cudakernel::check_error(const char * err_msg) {
#ifndef NDEBUG
  cudaDeviceSynchronize();
#endif
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    std::cout << "CUDA-ERROR " << err_msg << ": " << cudaGetErrorString(error) << std::endl;
  return (int) error;
}

