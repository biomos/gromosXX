#pragma once

#include <vector>
#include <stdexcept>

#include "gpu/cuda/cuheader.h"

namespace gpu {
    /**
     * @class CudaDeviceWorker
     * @brief Represents a single CUDA device and manages its resources.
     *
     * The CudaDeviceWorker class encapsulates the resources and operations associated
     * with a single CUDA device, such as streams, memory pools, and device-specific
     * kernel launches. It is designed to be used in a multi-GPU environment where
     * each GPU is managed independently.
     */
    class CudaDeviceWorker {
    public:
        /**
         * @brief Constructor for CudaDeviceWorker.
         * @param device_id The ID of the CUDA device to manage.
         */
        explicit CudaDeviceWorker(int device_id);

        /**
         * @brief Destructor for CudaDeviceWorker.
         *
         * Cleans up resources associated with the device, such as streams.
         */
        ~CudaDeviceWorker();

        /**
         * @brief Get the ID of the CUDA device managed by this worker.
         * @return The device ID.
         */
        int get_device_id() const;

        /**
         * @brief Get the CUDA stream associated with this device worker.
         * @return The CUDA stream.
         */
        CUSTREAM get_stream() const;

        /**
         * @brief Synchronize the device.
         *
         * Waits for all operations on the device to complete.
         */
        void synchronize() const;
#ifdef USE_CUDA
        /**
         * @brief Launch a kernel on this device.
         * @tparam KernelFunc The type of the kernel function.
         * @param kernel The kernel function to launch.
         * @param grid_dim The grid dimensions.
         * @param block_dim The block dimensions.
         * @param args The arguments to pass to the kernel.
         * @param shared_mem_size The amount of shared memory to allocate (default: 0).
         */
        template <typename KernelFunc, typename... Args>
        void launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size = 0);
#endif
    private:
        int device_id_; ///< The ID of the CUDA device managed by this worker.
        CUSTREAM stream_; ///< The CUDA stream associated with this device worker.
    };
}

// #include "cuda_device_worker.tcc" // Include template implementations