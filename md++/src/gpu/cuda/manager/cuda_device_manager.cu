
#include <sstream>

#include "gpu.h"

#include "gpu/cuda/cuheader.h"
#include "cuda_device_manager.h"

gpu::CudaDeviceManager::CudaDeviceManager() : device_count_(0), active_device_id_(-1) {
    CUERROR err = cudaGetDeviceCount(&device_count_);
    if (err != cudaSuccess || device_count_ == 0) {
        throw std::runtime_error("No CUDA devices available.");
    }
}

gpu::CudaDeviceManager::~CudaDeviceManager() {
    // No explicit cleanup needed for device manager
}

int gpu::CudaDeviceManager::get_device_count() const {
    return device_count_;
}

cudaDeviceProp gpu::CudaDeviceManager::get_device_properties(int device_id) const {
    if (device_id < 0 || device_id >= device_count_) {
        throw std::invalid_argument("Invalid device ID: " + std::to_string(device_id));
    }
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, device_id);
    return properties;
}

void gpu::CudaDeviceManager::set_device(int device_id) {
    if (device_id < 0 || device_id >= device_count_) {
        throw std::invalid_argument("Invalid device ID: " + std::to_string(device_id));
    }
    cudaSetDevice(device_id);
    active_device_id_ = device_id;
}

int gpu::CudaDeviceManager::get_active_device() const {
    return active_device_id_;
}

cudaStream_t gpu::CudaDeviceManager::create_stream() {
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    return stream;
}

void gpu::CudaDeviceManager::destroy_stream(cudaStream_t stream) {
    cudaStreamDestroy(stream);
}

std::string gpu::CudaDeviceManager::get_device_description(int device_id) const {
    if (device_id < 0 || device_id >= device_count_) {
        throw std::invalid_argument("Invalid device ID: " + std::to_string(device_id));
    }
    cudaDeviceProp properties = get_device_properties(device_id);

    std::ostringstream oss;
    oss << "Device " << device_id << ": " << properties.name << "\n"
        << "  Compute Capability: " << properties.major << "." << properties.minor << "\n"
        << "  Total Global Memory: " << (properties.totalGlobalMem / (1024 * 1024)) << " MB\n"
        << "  Multiprocessors: " << properties.multiProcessorCount << "\n"
        << "  Max Threads per Block: " << properties.maxThreadsPerBlock << "\n"
        << "  Max Threads per SM: " << properties.maxThreadsPerMultiProcessor << "\n";
    return oss.str();
}