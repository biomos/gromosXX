/**
 * @file reduction.h
 * @author capraz, poliak
 * device functions for efficient reduction
 */
#ifndef INCLUDED_REDUCTION_H
#define INCLUDED_REDUCTION_H

#include "types.h"

namespace gpu {
    /**
    * unroll last warp, no sync needed
    * continue summation when only the last warp is left (thread id <32). 
    * as instructions in one warp are synchronous no __syncthreads() needed.
    */
    template <unsigned int block_size, typename T>
    __device__  __forceinline__ void last_warp_reduce(volatile T *sdata, unsigned int tid);

    /**
    * Perform parallel reduction using sequential addressing of memory. 
    * First reads in multiple values into shared data of each block, then sums up all values within one block.
    */
    template <unsigned int block_size, typename T>
    __global__ void reduce(T *g_idata, T *g_odata, unsigned int g_size);

    /**
    * call kernel function reduce() iteratively until only one block is needed. 
    * device_input and device_output are switched after each iteration.
    * When only one block is needed iteration stops and the first element of device_output (the sum) is copied to host.
    */
    template <typename T, unsigned num_threads_per_block, unsigned num_elements_per_thread = 8>
    T calc_sum(unsigned int num_input_elements, T* device_input, T* device_output);

    /**
     * reduce variable in warp, Lane 0 gets the result
     */
    template <typename T>
    T reduce_in_warp(unsigned int num_input_elements, T* device_input, T* device_output);
}

#include "reduction.cu"

#endif