/**
 * @file gpu.cc
 * globals of the gpu library
 */

#include "../stdheader.h"
#include "config.h"

double gpu_ver = 0.10;

namespace gpu
{
  char const id[] = MD_VERSION;
  const char* get_id() { return id; }

#ifndef NDEBUG
  int debug_level = 0;
  int kernel_debug_level = 0;
  int cuda_debug_level = 0;
  int memory_debug_level = 0;
#endif
}