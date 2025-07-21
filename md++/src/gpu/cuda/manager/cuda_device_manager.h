
#pragma once

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <functional>

#include "cuda_memory_manager.h"
#include "cuda_device_worker.h"

namespace gpu {
    /**
     * @class CudaDeviceManager
     * @brief Manages a single CUDA device: initialization, memory, and kernel execution.
     * 
     * Responsible for:
     * - CUDA device context setup (cudaSetDevice)
     * - Memory allocation via CudaMemoryManager
     * - Task submission to device workers
     * - Synchronization and error checking
     */
    class CudaDeviceManager {
    public:
        /**
         * Constructor
         * @param device_id CUDA device ordinal (e.g. 0, 1, 2 ...)
         */
        explicit CudaDeviceManager(int device_id);

        /**
         * Destructor
         * Frees resources, synchronizes device
         */
        ~CudaDeviceManager();

        /**
         * @brief Disable copy
         */
        // CudaDeviceManager(const CudaDeviceManager&) = delete;
        // CudaDeviceManager& operator=(const CudaDeviceManager&) = delete;


        /**
         * @brief Allow shallow copy constructor, but warn
         */
        CudaDeviceManager(const CudaDeviceManager& other) {
            // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
            std::cerr << "Warning: Shallow copy of CudaDeviceManager at " << __FILE__
                    << ":" << __LINE__ << " in function " << __func__ << std::endl;
            this->m_device_name = other.m_device_name;
            this->m_workers = other.m_workers;
            this->m_device_id = other.m_device_id;
        }

        /**
         * @brief Allow shallow assignment operator, but warn
         */
        CudaDeviceManager& operator=(const CudaDeviceManager& other) {
            // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
            if (this != &other) {
                std::cerr << "Warning: Shallow copy assignment of CudaDeviceManager at " << __FILE__
                        << ":" << __LINE__ << " in function " << __func__ << std::endl;
                // Perform shallow copy
                this->m_device_name = other.m_device_name;
                this->m_workers = other.m_workers;
                this->m_device_id = other.m_device_id;
            }
            return *this;
        }

        /**
         * @brief Allow move
         */
        CudaDeviceManager(CudaDeviceManager&&) = default;
        CudaDeviceManager& operator=(CudaDeviceManager&&) = default;

        /**
         * Initialize the device context and resources
         * @return cudaError_t CUDA status
         */
        cudaError_t initialize();

        /**
         * @brief Get the properties of a specific CUDA device.
         * @param device_id The ID of the device.
         * @return The properties of the specified device.
         */
        cudaDeviceProp get_device_properties() const;

        /**
         * @brief Get a human-readable description of the CUDA device.
         * @return A string describing the device.
         */
        std::string get_device_description() const;


        /**
         * Allocates device memory of given size (in bytes)
         * @param size Number of bytes to allocate
         * @return void* pointer to device memory
         */
        void* allocate_memory(size_t size);

        /**
         * Frees previously allocated device memory
         * @param ptr Device pointer to free
         */
        void free_memory(void* ptr);

        /**
         * Submit a kernel task (or a callable that launches kernels) to this device.
         * Typically gets passed to a CudaDeviceWorker.
         * @param task Callable object (e.g. lambda) that executes CUDA kernels
         */
        void submit_task(std::function<void()> task);

        /**
         * @brief Synchronize the device (waits for all kernels to finish)
         * @throws std::runtime_error if CUDA fails.
         */
        void synchronize() const {
            CUDA_CHECK(cudaDeviceSynchronize());
        };

        
        /**
         * @brief Reset the active device to its default state.
         * @throws std::runtime_error if CUDA fails.
         */
        void reset() {
            CUDA_CHECK(cudaDeviceReset());
        };

        /**
         * Returns device ID
         */
        int get_device_id() const;

        /**
         * Get device name string
         */
        std::string device_name() const;

        /**
         * Adds a new worker
         */
        void add_worker();

        /**
         * Remove the last worker
         */
        void remove_worker();

        /**
         * Starts all worker threads
         */
        void start_all_workers(){
            for (auto& worker : m_workers) {
                worker.start();
            }
        }

        /**
         * Stops all worker threads gracefully
         */
        void stop_all_workers() {
            for (auto& worker : m_workers) {
                worker.stop();
            }
        }

    private:
        std::vector<CudaDeviceWorker> m_workers; // Map of device ID to workers
        std::mutex m_workers_mutex; // Protects access to next_worker
        std::string m_device_name;
        size_t m_next_worker;
        int m_device_id;
    };
}