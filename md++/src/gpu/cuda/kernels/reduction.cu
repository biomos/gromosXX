/**
 * @file reduction.cu
 * @author capraz, poliak
 * implementation device functions for efficient reduction
 */

#include "gpu.h"

#undef MODULE
#undef SUBMODULE
#define MODULE gpu
#define SUBMODULE math

#include <utility>

#include "reduction.h"

template <unsigned int block_size, typename T>
__device__  __forceinline__ void gpu::last_warp_reduce(volatile T *sdata, unsigned int tid) {
    if (block_size >=  64) sdata[tid] += sdata[tid + 32];
    if (block_size >=  32) sdata[tid] += sdata[tid + 16];
    if (block_size >=  16) sdata[tid] += sdata[tid +  8]; 
    if (block_size >=  8) sdata[tid] += sdata[tid +  4];
    if (block_size >=  4) sdata[tid] += sdata[tid +  2];
    if (block_size >=  2) sdata[tid] += sdata[tid +  1];
}

template <unsigned block_size, typename T>
__global__ void gpu::reduce(T *g_idata, T *g_odata, unsigned int g_size) {
    __shared__ T sdata[block_size];
    //T* volatile sdata = reinterpret_cast<T*>(shmem);
    unsigned tid = threadIdx.x;
    unsigned grid_size = block_size*2*gridDim.x;
    unsigned i = blockIdx.x*(block_size*2) + threadIdx.x;
    //initialize shared data with 0
    sdata[tid] = {0};

    //read global data into shared data, each thread sums up num_elements_per_thread while reading
    for (; i < g_size; i+=grid_size) {
        if (i+block_size < g_size){
            sdata[tid] += g_idata[i] + g_idata[i+block_size];
        } else {
            sdata[tid] += g_idata[i];
        }
    }
    __syncthreads();
    //sequential addressing within block with stride block_size/2; works up to NUM_THREADS_PER_BLOCK_REDUCTION == 1024
    if (block_size >= 512) {if (tid < 256) {sdata[tid] += sdata[tid + 256];} __syncthreads();}
    if (block_size >= 256) {if (tid < 128) {sdata[tid] += sdata[tid + 128];} __syncthreads();}
    if (block_size >= 128) {if (tid <  64) {sdata[tid] += sdata[tid +  64];} __syncthreads();}
    if (tid < 32) last_warp_reduce<block_size,T>(sdata, tid);
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <typename T, unsigned num_threads_per_block, unsigned num_elements_per_thread>
T gpu::calc_sum(unsigned int num_input_elements, T* device_input, T* device_output) {
    //NUM_THREADS_PER_BLOCK_REDUCTION has to be a power of 2 with maximum of 1024
    //number of elements each threads sums up when reading data from global to shared data
    //calculate the number of blocks needed depending on the block size, the size of the input and num_elements_per_thread
    unsigned num_blocks = (num_input_elements/num_threads_per_block)/num_elements_per_thread + 1;
    unsigned int size = num_input_elements;
    dim3 gridDim(num_blocks,1);
    dim3 blockDim(num_threads_per_block,1);

    while (true) {
        reduce<num_threads_per_block><<< gridDim, blockDim >>>(device_input, device_output, size);
        size = num_blocks;
        cudaDeviceSynchronize();
        if (num_blocks > 1){
            std::swap(device_input, device_output);
        } else {
            break;
        }
        num_blocks = (num_blocks/num_threads_per_block)/num_elements_per_thread+1;
    }
    //copy only the first element (i.e. the sum) of device_output to host_output
    T host_output;
    cudaMemcpy(&host_output, device_output, sizeof(T), cudaMemcpyDeviceToHost);
    return host_output;
}