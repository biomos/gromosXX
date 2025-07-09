

#include "cuda_device_worker.h"
#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"

gpu::CudaDeviceWorker::CudaDeviceWorker(int device_id) : device_id_(device_id), stream_(nullptr) {
    // Set the active device
    cudaSetDevice(device_id_);

    // Create a CUDA stream
    cudaStreamCreate(&stream_);
}

gpu::CudaDeviceWorker::~CudaDeviceWorker() {
    // Destroy the CUDA stream
    if (stream_) {
        cudaStreamDestroy(stream_);
    }
}

int gpu::CudaDeviceWorker::get_device_id() const {
    return device_id_;
}

cudaStream_t gpu::CudaDeviceWorker::get_stream() const {
    return stream_;
}

void gpu::CudaDeviceWorker::synchronize() const {
    cudaSetDevice(device_id_);
    cudaDeviceSynchronize();
}

template <typename KernelFunc, typename... Args>
void gpu::CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
    kernel<<<grid_dim, block_dim, shared_mem_size, stream_>>>(args...);

    // Check for errors after kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(err));
    }
}

// Explicit template instantiations (if needed)
template void gpu::CudaDeviceWorker::launch_kernel<void(*)(int*), int*>(void(*)(int*), dim3, dim3, int*, size_t);