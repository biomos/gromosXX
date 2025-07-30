/**
 * @file indexing.h
 * @author poliak
 * device function for easier indexing
 */

#pragma once

#define FULL_MASK 0xFFFFFFFF

__device__ __forceinline__ unsigned local_index() {
  return threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
}
__device__ __forceinline__ unsigned blocksize() {
  return blockDim.x * blockDim.y * blockDim.z;
}
__device__ __forceinline__ unsigned block_index() {
  return blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
}
__device__ __forceinline__ unsigned index() {
  return local_index() + blocksize() * block_index();
}

/**
 * fast integer square root
 */
__host__ __device__ unsigned fast_isqrt(const unsigned a) {
  unsigned g = 0x8000;
  unsigned c = 0x8000;
  for (;;) {
    if (g*g > a) {
      g ^= c;
    }
    c >>= 1;
    if (c == 0) {
      return g;
    }
    g |= c;
  }
}

