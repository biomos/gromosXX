

#include "cuda_device_worker.h"

template <typename KernelFunc, typename... Args>
void CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
    kernel<<<grid_dim, block_dim, shared_mem_size, stream_>>>(args...);

    // Check for errors after kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(err));
    }
}

// Explicit template instantiations (if needed)
template void CudaDeviceWorker::launch_kernel<void(*)(int*), int*>(void(*)(int*), dim3, dim3, int*, size_t);