#pragma once

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

/**
 * @class CudaMemoryManager
 * @brief Handles memory allocation, deallocation, and data transfers for CUDA.
 *
 * The CudaMemoryManager class provides a high-level interface for managing memory
 * on the GPU and transferring data between the host and device. It abstracts away
 * low-level CUDA memory management functions and ensures proper error handling.
 */
class CudaMemoryManager {
public:
    /**
     * @brief Allocate memory on the GPU.
     * @tparam T The type of the data to allocate.
     * @param size The number of elements to allocate.
     * @return A pointer to the allocated device memory.
     * @throws std::bad_alloc if the allocation fails.
     */
    template <typename T>
    T* allocate_device_memory(size_t size);

    /**
     * @brief Free memory on the GPU.
     * @tparam T The type of the data to free.
     * @param ptr A pointer to the device memory to free.
     */
    template <typename T>
    void free_device_memory(T* ptr);

    /**
     * @brief Allocate pinned (page-locked) memory on the host.
     * @tparam T The type of the data to allocate.
     * @param size The number of elements to allocate.
     * @return A pointer to the allocated host memory.
     * @throws std::bad_alloc if the allocation fails.
     */
    template <typename T>
    T* allocate_host_memory(size_t size);

    /**
     * @brief Free pinned (page-locked) memory on the host.
     * @tparam T The type of the data to free.
     * @param ptr A pointer to the host memory to free.
     */
    template <typename T>
    void free_host_memory(T* ptr);

    /**
     * @brief Copy data from the host to the device.
     * @tparam T The type of the data to copy.
     * @param device_ptr A pointer to the destination device memory.
     * @param host_ptr A pointer to the source host memory.
     * @param size The number of elements to copy.
     * @param stream The CUDA stream to use for the transfer (default: 0).
     */
    template <typename T>
    void copy_to_device(T* device_ptr, const T* host_ptr, size_t size, cudaStream_t stream = 0);

    /**
     * @brief Copy data from the device to the host.
     * @tparam T The type of the data to copy.
     * @param host_ptr A pointer to the destination host memory.
     * @param device_ptr A pointer to the source device memory.
     * @param size The number of elements to copy.
     * @param stream The CUDA stream to use for the transfer (default: 0).
     */
    template <typename T>
    void copy_to_host(T* host_ptr, const T* device_ptr, size_t size, cudaStream_t stream = 0);

    /**
     * @brief Copy data from one device memory location to another.
     * @tparam T The type of the data to copy.
     * @param dest_ptr A pointer to the destination device memory.
     * @param src_ptr A pointer to the source device memory.
     * @param size The number of elements to copy.
     * @param stream The CUDA stream to use for the transfer (default: 0).
     */
    template <typename T>
    void copy_device_to_device(T* dest_ptr, const T* src_ptr, size_t size, cudaStream_t stream = 0);
};

#include "cuda_memory_manager.tpp" // Include template implementations