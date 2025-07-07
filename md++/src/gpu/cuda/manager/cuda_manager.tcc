#pragma once

#include <cuda_runtime.h>
#include <stdexcept>

template <typename T>
gpu::cuvector<T> gpu::CudaManager::create_cuvector(int device_id, size_t size) {
    if (device_workers_.find(device_id) == device_workers_.end()) {
        throw std::runtime_error("Invalid device ID: " + std::to_string(device_id));
    }
    return gpu::cuvector<T>(size);
}

template <typename T>
void gpu::CudaManager::copy_to_device(int device_id, gpu::cuvector<T>& device_vector, const std::vector<T>& host_vector) {
    if (host_vector.size() > device_vector.size()) {
        throw std::runtime_error("Host vector size exceeds device vector size.");
    }
    memory_manager_.copy_to_device(device_vector.data(), host_vector.data(), host_vector.size(), get_stream(device_id));
}

template <typename T>
void gpu::CudaManager::copy_to_host(int device_id, std::vector<T>& host_vector, const gpu::cuvector<T>& device_vector) {
    if (host_vector.size() < device_vector.size()) {
        throw std::runtime_error("Host vector size is smaller than device vector size.");
    }
    memory_manager_.copy_to_host(host_vector.data(), device_vector.data(), device_vector.size(), get_stream(device_id));
}