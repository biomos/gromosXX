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

int cudakernel::checkError(const char * err_msg) {
#ifndef NDEBUG
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    std::cout << "CUDA-ERROR " << err_msg << ": " << cudaGetErrorString(error) << std::endl;
  return (int) error;
#else
  return 0;
#endif
}

