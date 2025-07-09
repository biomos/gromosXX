
#include <sstream>

#include "gpu/cuda/cuheader.h"
#include "cuda_device_manager.h"

gpu::CudaDeviceManager::CudaDeviceManager() : device_count_(0), active_device_id_(-1) {}

gpu::CudaDeviceManager::~CudaDeviceManager() {}

int gpu::CudaDeviceManager::get_device_count() const {
    return DISABLED(int);
}

CUDEVPROP gpu::CudaDeviceManager::get_device_properties(int device_id) const {
    return DISABLED(CUDEVPROP);
}

void gpu::CudaDeviceManager::set_device(int device_id) {
    DISABLED_VOID();
}

int gpu::CudaDeviceManager::get_active_device() const {
    return DISABLED(int);
}

CUSTREAM gpu::CudaDeviceManager::create_stream() {
    return DISABLED(CUSTREAM);
}

void gpu::CudaDeviceManager::destroy_stream(CUSTREAM stream) {
    DISABLED_VOID();
}

std::string gpu::CudaDeviceManager::get_device_description(int device_id) const {
    return DISABLED(std::string);
}