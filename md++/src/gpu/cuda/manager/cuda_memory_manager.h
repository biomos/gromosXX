#pragma once

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace gpu {
    /**
     * @class CudaMemoryManager
     * @brief Handles memory allocation, deallocation, and data transfers for CUDA devices.
     *
     * The CudaMemoryManager class provides a high-level interface for managing memory on CUDA devices.
     * It abstracts low-level CUDA memory operations and ensures efficient and error-resilient memory management.
     */
    class CudaMemoryManager {
    public:
        /**
         * @brief Constructor for CudaMemoryManager.
         */
        CudaMemoryManager();

        /**
         * @brief Destructor for CudaMemoryManager.
         */
        ~CudaMemoryManager();

        /**
         * @brief Allocate memory on the device.
         * @param size The size of the memory to allocate in bytes.
         * @return A pointer to the allocated device memory.
         * @throws std::runtime_error if memory allocation fails.
         */
        void* allocate_device_memory(size_t size);

        /**
         * @brief Free memory on the device.
         * @param device_ptr A pointer to the device memory to free.
         * @throws std::runtime_error if memory deallocation fails.
         */
        void free_device_memory(void* device_ptr);

        /**
         * @brief Allocate pinned host memory.
         * @param size The size of the memory to allocate in bytes.
         * @return A pointer to the allocated pinned host memory.
         * @throws std::runtime_error if memory allocation fails.
         */
        void* allocate_pinned_memory(size_t size);

        /**
         * @brief Free pinned host memory.
         * @param host_ptr A pointer to the pinned host memory to free.
         * @throws std::runtime_error if memory deallocation fails.
         */
        void free_pinned_memory(void* host_ptr);

        /**
         * @brief Copy data from host to device.
         * @param device_ptr A pointer to the destination device memory.
         * @param host_ptr A pointer to the source host memory.
         * @param size The size of the data to copy in bytes.
         * @throws std::runtime_error if the copy operation fails.
         */
        void copy_to_device(void* device_ptr, const void* host_ptr, size_t size);

        /**
         * @brief Copy data from device to host.
         * @param host_ptr A pointer to the destination host memory.
         * @param device_ptr A pointer to the source device memory.
         * @param size The size of the data to copy in bytes.
         * @throws std::runtime_error if the copy operation fails.
         */
        void copy_to_host(void* host_ptr, const void* device_ptr, size_t size);

        /**
         * @brief Copy data between two device memory locations.
         * @param dest_device_ptr A pointer to the destination device memory.
         * @param src_device_ptr A pointer to the source device memory.
         * @param size The size of the data to copy in bytes.
         * @throws std::runtime_error if the copy operation fails.
         */
        void copy_device_to_device(void* dest_device_ptr, const void* src_device_ptr, size_t size);

        /**
         * @brief Perform an asynchronous copy from host to device.
         * @param device_ptr A pointer to the destination device memory.
         * @param host_ptr A pointer to the source host memory.
         * @param size The size of the data to copy in bytes.
         * @param stream The CUDA stream to use for the asynchronous copy.
         * @throws std::runtime_error if the copy operation fails.
         */
        void async_copy_to_device(void* device_ptr, const void* host_ptr, size_t size, cudaStream_t stream);

        /**
         * @brief Perform an asynchronous copy from device to host.
         * @param host_ptr A pointer to the destination host memory.
         * @param device_ptr A pointer to the source device memory.
         * @param size The size of the data to copy in bytes.
         * @param stream The CUDA stream to use for the asynchronous copy.
         * @throws std::runtime_error if the copy operation fails.
         */
        void async_copy_to_host(void* host_ptr, const void* device_ptr, size_t size, cudaStream_t stream);

        /**
         * @brief Query the available and total memory on the device.
         * @param free_memory A reference to store the available memory in bytes.
         * @param total_memory A reference to store the total memory in bytes.
         * @throws std::runtime_error if the query operation fails.
         */
        void query_memory(size_t& free_memory, size_t& total_memory) const;

    private:
        /**
         * @brief Check the result of a CUDA operation and throw an exception if it failed.
         * @param result The result of the CUDA operation.
         * @param message A message to include in the exception if the operation failed.
         * @throws std::runtime_error if the CUDA operation failed.
         */
        void check_cuda_error(cudaError_t result, const std::string& message) const;
    };
}