#pragma once

#include <stdexcept>
#include <string>
#include <mutex>
#include <unordered_map>

#include "gpu/cuda/cuheader.h"

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
         * @brief Constructor
         */
        explicit CudaMemoryManager(int device_id);

        /**
         * @brief Destructor
         */
        ~CudaMemoryManager();

        /**
         * @brief Disable copy
         */
        CudaMemoryManager(const CudaMemoryManager&) = delete;
        CudaMemoryManager& operator=(const CudaMemoryManager&) = delete;

        /**
         * @brief Allow move
         */
        CudaMemoryManager(CudaMemoryManager&&) = default;
        CudaMemoryManager& operator=(CudaMemoryManager&&) = default;

        /**
         * @brief Allocate device memory; returns device pointer
         */
        void* allocate(std::size_t size_bytes);

        /**
         * @brief Free previously allocated device memory
         */
        void free(void* device_ptr);

        /**
         * @brief Initialize the memory manager.
         * @throws std::runtime_error.
         */
        void init();

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
        void copy_to_device(void* device_ptr, const void* host_ptr, size_t size, CUSTREAM stream);

        /**
         * @brief Copy data from device to host.
         * @param host_ptr A pointer to the destination host memory.
         * @param device_ptr A pointer to the source device memory.
         * @param size The size of the data to copy in bytes.
         * @throws std::runtime_error if the copy operation fails.
         */
        void copy_to_host(void* host_ptr, const void* device_ptr, size_t size, CUSTREAM stream);

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
        void async_copy_to_device(void* device_ptr, const void* host_ptr, size_t size, CUSTREAM stream);

        /**
         * @brief Perform an asynchronous copy from device to host.
         * @param host_ptr A pointer to the destination host memory.
         * @param device_ptr A pointer to the source device memory.
         * @param size The size of the data to copy in bytes.
         * @param stream The CUDA stream to use for the asynchronous copy.
         * @throws std::runtime_error if the copy operation fails.
         */
        void async_copy_to_host(void* host_ptr, const void* device_ptr, size_t size, CUSTREAM stream);

        
        // Synchronize the device to ensure all operations complete
        void synchronize();

        /**
         * @brief Query the available and total memory on the device.
         * @param free_memory A reference to store the available memory in bytes.
         * @param total_memory A reference to store the total memory in bytes.
         * @throws std::runtime_error if the query operation fails.
         */
        void query_memory(size_t& free_memory, size_t& total_memory) const;

        // Query currently allocated memory size
        std::size_t allocated_memory() const;

        // /**
        //  * @brief A specialized type, that allows user to create a transparent variable on the GPU
        //  * 
        //  * @tparam T Variable type
        //  */
        // template <typename T>
        // class Variable {
        // public:
        //     /**
        //      * @brief Constructor for Variable.
        //      * @param device_id The ID of the device to use.
        //      * @param file The source file where the Variable was created.
        //      * @param line The line number in the source file where the Variable was created.
        //      * @param func The function name where the Variable was created.
        //      */
        //     explicit Variable(int device_id = 0,
        //     const char* file = __FILE__,
        //     int line = __LINE__);
            
        //     ~Variable();

        //     Variable& operator=(const T& value);
        //     operator T() const;

        //     T* device_ptr();

        // private:
        //     T* m_device_data;

        //     void disabled(const char* file, int line, const char* func) const;
        // };


    private:
        int m_device_id;

        // Track allocations: pointer -> size
        std::unordered_map<void*, std::size_t> m_allocations;

        // Mutex for thread safety
        mutable std::mutex m_mutex;

        // lets play with vectors
        CUVECTOR_T<float> pos_umem;
        CUHVECTOR_T<float> pos_host;
        CUDVECTOR_T<float> pos_dev;
    };
}