#include "cuda_device_worker.h"

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