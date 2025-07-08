#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include "cuda_device_manager.h"
#include "cuda_device_worker.h"
#include "cuda_memory_manager.h"
#include "../memory/cuvector.h"

namespace gpu {
    /**
     * @class CudaManager
     * @brief High-level orchestrator for multi-GPU CUDA operations in a molecular dynamics program.
     *
     * The CudaManager class initializes the CUDA environment, manages multiple devices, and coordinates
     * memory and kernel execution across GPUs. It provides a clean interface for interacting with CUDA
     * resources and abstracts away low-level details.
     */
    class CudaManager {
    public:
        /**
         * @brief Constructor for CudaManager.
         */
        CudaManager();

        /**
         * @brief Destructor for CudaManager.
         */
        ~CudaManager();

        /**
         * @brief Initialize the CUDA environment and select devices.
         * @param device_ids A vector of device IDs to use. If empty, all available devices are used.
         * @throws std::runtime_error if no devices are available or initialization fails.
         */
        void initialize(const std::vector<int>& device_ids = {});

        /**
         * @brief Get the number of active GPUs.
         * @return The number of active GPUs.
         */
        size_t get_device_count() const;

        /**
         * @brief Get the CUDA stream for a specific device.
         * @param device_id The ID of the device.
         * @return The CUDA stream for the specified device.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        cudaStream_t get_stream(int device_id) const;

        /**
         * @brief Synchronize all devices.
         */
        void synchronize_all();

        /**
         * @brief Synchronize a specific device.
         * @param device_id The ID of the device to synchronize.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        void synchronize_device(int device_id);

        /**
         * @brief Create a cuvector for managing device memory on a specific GPU.
         * @tparam T The type of elements in the cuvector.
         * @param device_id The ID of the device.
         * @param size The number of elements to allocate.
         * @return A cuvector of the specified size.
         * @throws std::invalid_argument if the device ID is invalid.
         * @throws std::runtime_error if memory allocation fails.
         */
        template <typename T>
        cuvector<T> create_cuvector(int device_id, size_t size);

        /**
         * @brief Copy data from a host vector to a cuvector on a specific GPU.
         * @tparam T The type of elements in the vectors.
         * @param device_id The ID of the device.
         * @param device_vector The cuvector on the device.
         * @param host_vector The host vector containing the data.
         * @throws std::invalid_argument if the device ID is invalid.
         * @throws std::runtime_error if the copy operation fails.
         */
        template <typename T>
        void copy_to_device(int device_id, cuvector<T>& device_vector, const std::vector<T>& host_vector);

        /**
         * @brief Copy data from a cuvector on a specific GPU to a host vector.
         * @tparam T The type of elements in the vectors.
         * @param device_id The ID of the device.
         * @param host_vector The host vector to receive the data.
         * @param device_vector The cuvector on the device.
         * @throws std::invalid_argument if the device ID is invalid.
         * @throws std::runtime_error if the copy operation fails.
         */
        template <typename T>
        void copy_to_host(int device_id, std::vector<T>& host_vector, const cuvector<T>& device_vector);

        /**
         * @brief Get a human-readable description of all active devices.
         * @return A vector of strings describing the active devices.
         */
        std::vector<std::string> get_active_device_descriptions() const;

    private:
        /**
         * @brief Validate a device ID.
         * @param device_id The ID of the device to validate.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        void validate_device_id(int device_id) const;

        std::unique_ptr<CudaDeviceManager> device_manager_; ///< Manages CUDA devices.
        std::unordered_map<int, std::unique_ptr<CudaDeviceWorker>> device_workers_; ///< Workers for each active device.
        CudaMemoryManager memory_manager_; ///< Handles memory allocation and transfers.
    };
}

#include "cuda_manager.tcc" // Include template implementations