
#include "gpu.h"

#include "cuda_device_worker.h"
#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"

gpu::CudaDeviceWorker::CudaDeviceWorker(int device_id) : m_stream(nullptr), m_device_id(device_id) {
    // Set the active device
    CUDA_CHECK(cudaSetDevice(m_device_id));

    // Create a CUDA stream
    CUDA_CHECK(cudaStreamCreate(&m_stream));
}

gpu::CudaDeviceWorker::~CudaDeviceWorker() {
    // Destroy the CUDA stream
    if (m_stream) {
        CUDA_CHECK(cudaStreamDestroy(m_stream));
    }
}

int gpu::CudaDeviceWorker::get_device_id() const {
    return m_device_id;
}

cudaStream_t gpu::CudaDeviceWorker::get_stream() const {
    return m_stream;
}

void gpu::CudaDeviceWorker::synchronize() const {
    CUDA_CHECK(cudaDeviceSynchronize());
}

template <typename KernelFunc, typename... Args>
void gpu::CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
    kernel<<<grid_dim, block_dim, shared_mem_size, m_stream>>>(args...);

    // Check for errors after kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(err));
    }
}

// Explicit template instantiations (if needed)
// template void gpu::CudaDeviceWorker::launch_kernel<void(*)(int*), int*>(void(*)(int*), dim3, dim3, int*, size_t);

// template <typename KernelFunc, typename... Args>
// void gpu::CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
//     // Set the active device
//     cudaSetDevice(device_id_);

//     // Launch the kernel
//     kernel<<<grid_dim, block_dim, shared_mem_size, m_stream>>>(args...);

//     // Check for kernel launch errors
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess) {
//         throw std::runtime_error("Kernel launch failed: " + std::string(cudaGetErrorString(err)));
//     }
// }