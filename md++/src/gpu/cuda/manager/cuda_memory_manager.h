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
         * @brief Allocate raw device memory; returns device pointer
         */
        void* allocate(std::size_t size_bytes);

        /**
         * @brief Allocate device memory for custom types
         */
        template <typename T>
        T* allocate(const T& host_data) {
            cudaSetDevice(m_device_id);
            std::lock_guard<std::mutex> lock(m_mutex);
            T* devptr;
            cudaMalloc(&devptr, sizeof(T));
            cudaMemcpy(devptr, &host_data, sizeof(T), cudaMemcpyHostToDevice);
            m_allocations[devptr] = sizeof(T);
            return devptr;
        };

        /**
         * @brief Allocate unified memory for custom types
         */
        template <typename T>
        T* allocate_managed(const T& host_data) {
            std::lock_guard<std::mutex> lock(m_mutex);
            T* devptr;
            cudaMallocManaged(&devptr, sizeof(T));
            *devptr = host_data;
            m_allocations[devptr] = sizeof(T);
            return devptr;
        };

        /**
         * @brief Allocate device memory for custom arrays
         */
        template <typename T>
        T* allocate_array(size_t count) {
            std::lock_guard<std::mutex> lock(m_mutex);

            T* devptr;
            cudaMalloc(&devptr, count * sizeof(T));
            m_allocations[devptr] = count * sizeof(T);

            return devptr;
        };

        // Deallocate memory and remove it from the tracking map
        template <typename T>
        void deallocate(T* device_ptr) {
            std::lock_guard<std::mutex> lock(m_mutex);

            auto it = m_allocations.find(device_ptr);
            if (it == m_allocations.end()) {
                throw std::runtime_error("Attempt to free memory that was not allocated by CudaMemoryManager.");
            }

            cudaFree(device_ptr);
            m_allocations.erase(it);
        }

        // Copy data to device memory
        template <typename T>
        void copyToDevice(T* devicePtr, const T& hostData) {
            std::lock_guard<std::mutex> lock(m_mutex);

            auto it = m_allocations.find(devicePtr);
            if (it == m_allocations.end()) {
                throw std::runtime_error("Attempt to copy to memory that was not allocated by CudaMemoryManager.");
            }

            if (it->second < sizeof(T)) {
                throw std::runtime_error("Insufficient memory allocated for the copy operation.");
            }

            cudaMemcpy(devicePtr, &hostData, sizeof(T), cudaMemcpyHostToDevice);
        }

        // Copy data from device memory
        template <typename T>
        void copyToHost(const T* devicePtr, T& hostData) {
            std::lock_guard<std::mutex> lock(m_mutex);

            auto it = m_allocations.find(const_cast<T*>(devicePtr));
            if (it == m_allocations.end()) {
                throw std::runtime_error("Attempt to copy from memory that was not allocated by CudaMemoryManager.");
            }

            if (it->second < sizeof(T)) {
                throw std::runtime_error("Insufficient memory allocated for the copy operation.");
            }

            cudaMemcpy(&hostData, devicePtr, sizeof(T), cudaMemcpyDeviceToHost);
        }

        // Debugging: Print all tracked allocations
        void printAllocations() const {
            std::lock_guard<std::mutex> lock(m_mutex);

            for (const auto& [ptr, size] : m_allocations) {
                std::cout << "Device pointer: " << ptr << ", Size: " << size << " bytes\n";
            }
        }

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

        // Debugging: Print all tracked allocations
        void print_allocations() const {
            std::lock_guard<std::mutex> lock(m_mutex);

            for (const auto& [ptr, size] : m_allocations) {
                std::cout << "Device pointer: " << ptr << ", Size: " << size << " bytes\n";
            }
        }

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
    };
}