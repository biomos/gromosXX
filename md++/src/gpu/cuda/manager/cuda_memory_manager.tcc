#pragma once

#include <cuda_runtime.h>
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
    cudaError_t err = cudaFreeHost(ptr);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to free host memory: " + std::string(cudaGetErrorString(err)));
    }
}

template <typename T>
void CudaMemoryManager::copy_to_device(T* device_ptr, const T* host_ptr, size_t size, cudaStream_t stream) {
    cudaError_t err = cudaMemcpyAsync(device_ptr, host_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to copy data to device: " + std::string(cudaGetErrorString(err)));
    }
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