/**
 * @file gpu.h
 * globals of the gpu library
 */

#pragma once

extern double gpu_ver;

namespace gpu {
  extern const char id[];
  const char* get_id();

#ifndef NDEBUG
  extern int debug_level;
  extern int kernel_debug_level;
  extern int cuda_debug_level;
  extern int memory_debug_level;
#endif

} // namespace gpu
