#pragma once

#include <stdexcept>

#include "gpu/cuda/cuheader.h"


/*
    Possible cuvector<T> implementation
*/
// template <template<typename, typename> class VecT, typename T, typename Alloc /* = gpu::CuMAllocator<T> */>
// VecT<T, Alloc> gpu::CudaManager::create_cuvector(int device_id, size_t size) {
//     // static_assert(std::is_same_v<VecT<T, Alloc>, gpu::cuvector<T>> ||
//     //               std::is_same_v<VecT<T, Alloc>, gpu::cuhvector<T>> ||
//     //               std::is_same_v<VecT<T, Alloc>, gpu::cudvector<T>>,
//     //               "Unsupported vector type");

//     if (device_workers_.find(device_id) == device_workers_.end()) {
//         throw std::runtime_error("Invalid device ID");
//     }

//     return VecT<T, Alloc>(size);
// }

template <typename T>
gpu::CUVECTOR_T<T> gpu::CudaManager::create_cuvector(int device_id, size_t size) {
    if (device_workers_.find(device_id) == device_workers_.end()) {
        throw std::runtime_error("Invalid device ID: " + std::to_string(device_id));
    }
    return gpu::CUVECTOR_T<T>(size);
}

template <typename T>
void gpu::CudaManager::copy_to_device(int device_id, gpu::CUVECTOR_T<T>& device_vector, const std::vector<T>& host_vector) {
    if (host_vector.size() > device_vector.size()) {
        throw std::runtime_error("Host vector size exceeds device vector size.");
    }
    memory_manager_.copy_to_device(device_vector.data(), host_vector.data(), host_vector.size(), get_stream(device_id));
}

template <typename T>
void gpu::CudaManager::copy_to_host(int device_id, std::vector<T>& host_vector, const gpu::CUVECTOR_T<T>& device_vector) {
    if (host_vector.size() < device_vector.size()) {
        throw std::runtime_error("Host vector size is smaller than device vector size.");
    }
    memory_manager_.copy_to_host(host_vector.data(), device_vector.data(), device_vector.size(), get_stream(device_id));
}

/**
 * @brief Construct a new gpu::Cuda Manager::Variable<T>::Variable object
 * 
 * @tparam T 
 * @param device_id 
 * @param file 
 * @param line 
 * @param func 
 */
template <typename T>
gpu::CudaManager::Variable<T>::Variable(int device_id, const char* file, int line)
    : device_id_(device_id), device_data_(nullptr) {
#ifdef USE_CUDA
    cudaSetDevice(device_id_);
    cudaMalloc(&device_data_, sizeof(T));
#else
    CUDA_VARIABLE_DISABLED();
#endif
}

template <typename T>
gpu::CudaManager::Variable<T>::~Variable() {
#ifdef USE_CUDA
    cudaFree(device_data_);
#endif
}

template <typename T>
typename gpu::CudaManager::Variable<T>& gpu::CudaManager::Variable<T>::operator=(const T& value) {
#ifdef USE_CUDA
    cudaMemcpy(device_data_, &value, sizeof(T), cudaMemcpyHostToDevice);
    return *this;
#else
    CUDA_VARIABLE_DISABLED();
    return *this;
#endif
}

template <typename T>
gpu::CudaManager::Variable<T>::operator T() const {
    T host_copy{};
#ifdef USE_CUDA
    cudaMemcpy(&host_copy, device_data_, sizeof(T), cudaMemcpyDeviceToHost);
    return host_copy;
#else
    CUDA_VARIABLE_DISABLED();
    return host_copy;
#endif
}

template <typename T>
T* gpu::CudaManager::Variable<T>::device_ptr() {
#ifdef USE_CUDA
    return device_data_;
#else
    CUDA_VARIABLE_DISABLED();
    return nullptr;
#endif
}

template <typename T>
void gpu::CudaManager::Variable<T>::disabled(const char* file, int line, const char* func) const {
    std::ostringstream oss;
    oss << "CUDA is disabled in this build.\n"
        << "Attempted to use gpu::CudaManager::variable<" << typeid(T).name() << "> at:\n"
        << "  " << file << ":" << line << " (" << func << ")";
    throw std::runtime_error(oss.str());
}
