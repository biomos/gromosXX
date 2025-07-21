#pragma once

#ifdef USE_CUDA
    #include <cuda_runtime.h>
#endif
#include <stdexcept>
#include <string>

template <typename T>
T* CudaMemoryManager::allocate_device_memory(size_t size) {
    T* ptr = nullptr;
    cudaError_t err = cudaMalloc(&ptr, size * sizeof(T));
    if (err != cudaSuccess) {
        throw std::bad_alloc();
    }
    return ptr;
}

template <typename T>
void CudaMemoryManager::free_device_memory(T* ptr) {
    cudaError_t err = cudaFree(ptr);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to free device memory: " + std::string(cudaGetErrorString(err)));
    }
}

template <typename T>
T* CudaMemoryManager::allocate_host_memory(size_t size) {
    T* ptr = nullptr;
    cudaError_t err = cudaHostAlloc(&ptr, size * sizeof(T), cudaHostAllocDefault);
    if (err != cudaSuccess) {
        throw std::bad_alloc();
    }
    return ptr;
}

template <typename T>
void CudaMemoryManager::free_host_memory(T* ptr) {
    CUDA_CHECK(cudaFreeHost(ptr));
}

template <typename T>
void CudaMemoryManager::copy_to_device(T* device_ptr, const T* host_ptr, size_t size, cudaStream_t stream) {
    CUDA_CHECK(cudaMemcpyAsync(device_ptr, host_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream));
}

template <typename T>
void CudaMemoryManager::copy_to_host(T* host_ptr, const T* device_ptr, size_t size, cudaStream_t stream) {
    cudaError_t err = cudaMemcpyAsync(host_ptr, device_ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy data to host: " + std::string(cudaGetErrorString(err)));
    }
}

template <typename T>
void CudaMemoryManager::copy_device_to_device(T* dest_ptr, const T* src_ptr, size_t size, cudaStream_t stream) {
    cudaError_t err = cudaMemcpyAsync(dest_ptr, src_ptr, size * sizeof(T), cudaMemcpyDeviceToDevice, stream);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy data between device memory: " + std::string(cudaGetErrorString(err)));
    }
}

/**
 * class Variable template definitions
 */


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