#pragma once

#include <cuda_runtime.h>
#include <stdexcept>

template <typename KernelFunc, typename... Args>
void CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
    // Set the active device
    cudaSetDevice(device_id_);

    // Launch the kernel
    kernel<<<grid_dim, block_dim, shared_mem_size, stream_>>>(args...);

    // Check for kernel launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error("Kernel launch failed: " + std::string(cudaGetErrorString(err)));
    }
}