#include "cuda_device_worker.h"
#include "gpu/cuda/cuheader.h"

gpu::CudaDeviceWorker::CudaDeviceWorker(int device_id) : device_id_(device_id), stream_(nullptr) {}

gpu::CudaDeviceWorker::~CudaDeviceWorker() {}

int gpu::CudaDeviceWorker::get_device_id() const {
    return DISABLED(int);
}

CUSTREAM gpu::CudaDeviceWorker::get_stream() const {
    return DISABLED(CUSTREAM);
}

void gpu::CudaDeviceWorker::synchronize() const {
    DISABLED_VOID();
}