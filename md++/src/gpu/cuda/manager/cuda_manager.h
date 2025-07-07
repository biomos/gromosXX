#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include "cuda_device_manager.h"
#include "cuda_device_worker.h"
#include "cuda_memory_manager.h"
#include "cuvector.h"

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
     */
    cudaStream_t get_stream(int device_id) const;

    /**
     * @brief Create a cuvector for managing device memory on a specific GPU.
     * @tparam T The type of elements in the cuvector.
     * @param device_id The ID of the device.
     * @param size The number of elements to allocate.
     * @return A cuvector of the specified size.
     */
    template <typename T>
    cuvector<T> create_cuvector(int device_id, size_t size);

    /**
     * @brief Copy data from a host vector to a cuvector on a specific GPU.
     * @tparam T The type of elements in the vectors.
     * @param device_id The ID of the device.
     * @param device_vector The cuvector on the device.
     * @param host_vector The host vector containing the data.
     */
    template <typename T>
    void copy_to_device(int device_id, cuvector<T>& device_vector, const std::vector<T>& host_vector);

    /**
     * @brief Copy data from a cuvector on a specific GPU to a host vector.
     * @tparam T The type of elements in the vectors.
     * @param device_id The ID of the device.
     * @param host_vector The host vector to receive the data.
     * @param device_vector The cuvector on the device.
     */
    template <typename T>
    void copy_to_host(int device_id, std::vector<T>& host_vector, const cuvector<T>& device_vector);

    /**
     * @brief Synchronize all devices.
     */
    void synchronize_all();

private:
    std::unique_ptr<CudaDeviceManager> device_manager_; ///< Manages CUDA devices.
    std::unordered_map<int, std::unique_ptr<CudaDeviceWorker>> device_workers_; ///< Workers for each active device.
    CudaMemoryManager memory_manager_; ///< Handles memory allocation and transfers.
};

#include "cuda_manager.tpp" // Include template implementations