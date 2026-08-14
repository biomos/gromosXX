#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include "gpu/cuda/cuheader.h"

#ifdef USE_CUDA
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#endif

namespace topology {
    class Topology;
}
namespace configuration {
    class Configuration;
}

#define CUDA_VARIABLE_DISABLED() disabled(__FILE__, __LINE__, __func__)
namespace gpu {
    class CudaDeviceManager;
    // class CudaMemoryManager;
    // class CudaDeviceWorker;
    /**
     * @class CudaManager
     * @brief High-level orchestrator for multi-GPU CUDA operations in a molecular dynamics program.
     *
     * The CudaManager class initializes the CUDA environment, manages multiple devices, and coordinates
     * memory and kernel execution across GPUs. It provides a clean interface for interacting with CUDA
     * resources and abstracts away low-level details.
     * 
     * This is the only interface for host code to access CUDA, to encapsulate all cuda code separately
     * and not expose it to the CPU-only code.
     * Function calls, that submit to GPU should go only and only over this manager.
     * Probably should be called CudaInterface, we will see.
     * 
     */

    #define CUDA_MANAGER_COPY_WARNING \
    _Pragma("message(\"Warning: Shallow copy of CudaManager detected.\")")

    class CudaManager {
        public:
            /**
             * @brief Constructor
             */
            CudaManager();

            /**
             * Destructor
             */
            // ~CudaManager();

            /**
             * @brief Disable copy
             */
            // CudaManager(const CudaManager&) = delete;
            // CudaManager& operator=(const CudaManager&) = delete;

            /**
             * @brief Allow shallow copy constructor, but warn
             */
            CudaManager(const CudaManager& other);

            /**
             * @brief Allow shallow assignment operator, but warn
             */
            CudaManager& operator=(const CudaManager& other);


            /**
             * @brief Allow move
             */
            CudaManager(CudaManager&&) = default;
            CudaManager& operator=(CudaManager&&) = default;

            /**
             * @brief Initialize the CUDA environment and select devices.
             * @param device_ids A vector of device IDs to use. If empty, all available devices are used.
             * @throws std::runtime_error if no devices are available or initialization fails.
             */
            void init(const std::vector<int>& device_ids = {});

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
            // CUSTREAM get_stream(int device_id) const;

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

            // /**
            //  * @brief Allocate raw device memory; returns device pointer
            //  */
            // void* allocate(std::size_t size_bytes);

            // /**
            //  * @brief Allocate device memory for custom types
            //  */
            // template <typename T>
            // T* allocate(const T& host_data) {
            //     for (const auto& [device_id, device_manager] : m_device_managers) {
            //         dm->memory().allocate(host_data);
            //     }
            //     std::lock_guard<std::mutex> lock(m_mutex);
            //     T* devptr;
            //     cudaMalloc(&devptr, sizeof(T));
            //     cudaMemcpy(devptr, &host_data, sizeof(T), cudaMemcpyHostToDevice);
            //     m_allocations[devptr] = sizeof(T);
            //     return devptr;
            // };

            /**
             * @brief Create a cuvector for managing device memory on a specific GPU.
             * @tparam T The type of elements in the cuvector.
             * @param device_id The ID of the device.
             * @param size The number of elements to allocate.
             * @return A cuvector of the specified size.
             * @throws std::invalid_argument if the device ID is invalid.
             * @throws std::runtime_error if memory allocation fails.
             */
            // template <template<typename, typename> class VecT, typename T, typename Alloc /* = gpu::CuMAllocator<T> */>
            // VecT<T, Alloc> create_cuvector(int device_id, size_t size);
            // template <typename T>
            // gpu::CUVECTOR_T<T> create_cuvector(int device_id, size_t size);

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
            void copy_to_device(int device_id, gpu::CUVECTOR_T<T>& device_vector, const std::vector<T>& host_vector);

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
            void copy_to_host(int device_id, std::vector<T>& host_vector, const gpu::CUVECTOR_T<T>& device_vector);

            /**
             * @brief Get a human-readable description of all active devices.
             * @return A vector of strings describing the active devices.
             */
            std::vector<std::string> get_active_device_descriptions() const;

            /**
             * @brief Automatically select the best CUDA device based on properties.
             * @return The ID of the selected device.
             * @throws std::runtime_error if no suitable device is found.
             */
            int select_best_device() const;

#ifdef USE_CUDA
            /**
             * @brief Identity-keyed GPU mirror cache for topology::Topology
             * (PLAN.md §3.2). Builds the mirror on first call for a given
             * topo.id(); returns the cached one on subsequent calls unless
             * `force_resync` is set -- topology data is static for a
             * normal run, so this is rare (lambda/perturbation-topology
             * updates are the exception). A cache entry that belonged to a
             * different, no-longer-live Topology which happened to be
             * destroyed and have its heap address reused is never a risk
             * here: the cache is keyed on `id()`, a process-wide token
             * that's never reused, not on the object's address -- an
             * unrelated object at the same address has a different id and
             * is correctly treated as a cache miss, not a stale hit.
             */
            gpu::Topology::View topology_view(const topology::Topology & topo,
                                               bool force_resync = false);

            /**
             * @brief Identity-keyed GPU mirror cache for
             * configuration::Configuration (PLAN.md §3.2). Builds the
             * mirror (full sync -- pos/vel/force/constraint_force, current
             * and old, plus box/tensors) on first call for a given
             * conf.id(). Subsequent calls default to a cheap
             * positions+velocities-only resync (`sync_pos_vel = true`)
             * since those change every step, unlike the rest of the
             * mirrored state; pass `false` if the caller already knows
             * positions/velocities haven't changed since the last call
             * (e.g. a second read within the same step).
             */
            gpu::Configuration::View configuration_view(configuration::Configuration & conf,
                                                          bool sync_pos_vel = true,
                                                          bool full_resync = false);

            /**
             * @brief Publish the GPU mirror's current positions/velocities
             * back to the CPU-authoritative configuration::Configuration.
             * For GPU-native integrators (Leap_Frog_*<gpuBackend>) that
             * leave their result resident on the GPU mirror across
             * multiple algorithms and only need one sync-back at the very
             * end, rather than after every kernel. No-op if `conf` has
             * never been mirrored (nothing to publish).
             */
            void sync_configuration_from_device(configuration::Configuration & conf);
#endif

        private:
            /**
             * @brief Validate a device ID.
             * @param device_id The ID of the device to validate.
             * @throws std::invalid_argument if the device ID is invalid.
             */
            void validate_device_id(int device_id) const;
#ifdef USE_CUDA
            std::unordered_map<int, std::shared_ptr<CudaDeviceManager> > m_device_managers; ///< Managers for each active device.

            std::unordered_map<std::size_t, std::unique_ptr<gpu::Topology> > m_topologies;
            std::unordered_map<std::size_t, std::unique_ptr<gpu::Configuration> > m_configurations;

            // 1-entry fast path for the overwhelmingly common single-
            // topology/single-configuration case, avoiding a hashmap
            // lookup on every call. 0 is never a real id (util::
            // next_identity_token() starts at 1), so it's a safe "nothing
            // cached yet" sentinel.
            std::size_t m_last_topo_id = 0;
            gpu::Topology * m_last_topo_gpu = nullptr;
            std::size_t m_last_conf_id = 0;
            gpu::Configuration * m_last_conf_gpu = nullptr;
#endif
    };
}
