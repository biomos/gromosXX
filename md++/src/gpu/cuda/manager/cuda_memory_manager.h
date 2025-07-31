#pragma once

#include <stdexcept>
#include <string>
#include <mutex>
#include <unordered_map>

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/cupointer.h"

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

        /** Allow allocating only with smart pointers
         * we track all allocations and allow deallocation only of the pointers we own
            1. allocate<T>(size_t N) - allocates array of type T of size N
            2. copy_to_device<T>(T& var) - allocate on device and copy (only for primitive types)
            3. copy_to_device<T>(vector_type<T>& vec) - allocate on device and copy (only for vector types)
            4. copy_to_device<T>(T& source, cu_unique_ptr<T> dest_pointer) - copy to existing pointer (only for primitive types)
            5. copy_to_device<T>(vector_type<T>& source_vec, cu_unique_ptr<T> dest_pointer) - copy to existing pointer (only for vector types)
            6. // copy_to_host variants
            7. deallocate<T>(cu_unique_ptr<T> T)
        */

        template <typename T>
        using cumm_unique_ptr = std::unique_ptr<T, CuDeleter<T, CudaMemoryManager>>;

        /**
         * @brief Allocate raw device memory; returns device pointer
         * We try not allowing raw pointers
         */
        // template <typename T>
        // T* allocate(std::size_t size_bytes);

        /**
         * @brief Allocate smart device memory for custom types
         */
        template <typename T>
        cumm_unique_ptr<T> allocate(size_t count = 1);

        /**
         * @brief allocate and copy to device memory for custom types
         */
        template<typename T>
        cumm_unique_ptr<std::remove_cv_t<T>> copy_to_device(const T& src);

        /**
         * @brief allocate and copy to device memory for containers
         * 
         * @tparam C source vector type
         * @tparam T source vector value_type
         * @tparam U destination type (for e.g. casting double to float)
         * @param vec 
         * @return cumm_unique_ptr<T> 
         */
        template<typename C,
                typename T = typename C::value_type,
                typename U = T>
        cumm_unique_ptr<U> copy_to_device(const C& vec);

        /**
         * @brief Allocate unified memory for custom types
         */
        template <typename T>
        T* allocate_managed(const T& host_data);

        // /**
        //  * @brief Allocate device memory for custom arrays
        //  */
        // template <typename T>
        // T* allocate_array(size_t count) {
        //     std::lock_guard<std::mutex> lock(m_mutex);

        //     T* devptr;
        //     cudaMalloc(&devptr, count * sizeof(T));
        //     m_allocations[devptr] = count * sizeof(T);

        //     return devptr;
        // };

        // Deallocate memory and remove it from the tracking map - only for raw pointers
        void deallocate(void* device_ptr);

        void validate(void* ptr) const;

        // Debugging: Print all tracked allocations
        void print_allocations() const;

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
        std::unordered_map<void*, size_t> m_allocations;

        // Mutex for thread safety
        mutable std::mutex m_mutex;
    };
}