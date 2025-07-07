#include "cuda_device_worker.h"

CudaDeviceWorker::CudaDeviceWorker(int device_id) : device_id_(device_id), stream_(nullptr) {
    // Set the active device
    cudaSetDevice(device_id_);

    // Create a CUDA stream
    cudaStreamCreate(&stream_);
}

CudaDeviceWorker::~CudaDeviceWorker() {
    // Destroy the CUDA stream
    if (stream_) {
        cudaStreamDestroy(stream_);
    }
}

int CudaDeviceWorker::get_device_id() const {
    return device_id_;
}

cudaStream_t CudaDeviceWorker::get_stream() const {
    return stream_;
}

void CudaDeviceWorker::synchronize() const {
    cudaSetDevice(device_id_);
    cudaDeviceSynchronize();
}