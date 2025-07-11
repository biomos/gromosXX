/**
 * @file gpu.cc
 * globals of the gpu library
 */

#include "../stdheader.h"
#include "config.h"
#include "gpu.h"

double gpu_ver = 0.10;

namespace gpu {
  const char id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int cuda_debug_level = 0; // CUDA manager debug
  int device_debug_level = 0; // CUDA device manager debug 
  int memory_debug_level = 0; // CUDA memory manager debug 
  int kernel_debug_level = 0; // CUDA kernel debug 
#endif
}
