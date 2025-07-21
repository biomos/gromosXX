#pragma once

// #include <vector>
// #include <stdexcept>
#include <queue>
#include <condition_variable>
#include <functional>
#include <thread>

#include "gpu/cuda/cuheader.h"

namespace gpu {
    /**
     * @class CudaDeviceWorker
     * @brief Represents a single CUDA stream to perform tasks on demand.
     *
     * The CudaDeviceWorker class encapsulates the stream, task queue and submission
     * within a single CUDA device.
     */

    class CudaDeviceWorker {
    public:
        /**
         * Constructor
         * @param device_id The CUDA device ID this worker is associated with.
         */
        explicit CudaDeviceWorker(int device_id);

        /**
         * Destructor
         * Ensures the worker thread is stopped and resources are cleaned up.
         */
        ~CudaDeviceWorker();

        /**
         * @brief Disable copy
         */
        // CudaDeviceWorker(const CudaDeviceWorker&) = delete;
        // CudaDeviceWorker& operator=(const CudaDeviceWorker&) = delete;

        /**
         * @brief Allow shallow copy constructor, but warn
         */
        CudaDeviceWorker(const CudaDeviceWorker& other) {
            // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
            std::cerr << "Warning: Shallow copy of CudaDeviceWorker at " << __FILE__
                    << ":" << __LINE__ << " in function " << __func__ << std::endl;
            this->m_task_queue = other.m_task_queue;
            this->m_stream = other.m_stream;
            // this->m_task_queue = other.m_task_queue;
            this->m_device_id = other.m_device_id;
            this->m_stop_worker = other.m_stop_worker;
        }

        /**
         * @brief Allow shallow assignment operator, but warn
         */
        CudaDeviceWorker& operator=(const CudaDeviceWorker& other) {
            // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
            if (this != &other) {
                std::cerr << "Warning: Shallow copy assignment of CudaDeviceWorker at " << __FILE__
                        << ":" << __LINE__ << " in function " << __func__ << std::endl;
                // Perform shallow copy
                this->m_task_queue = other.m_task_queue;
                this->m_stream = other.m_stream;
                // this->m_task_queue = other.m_task_queue;
                this->m_device_id = other.m_device_id;
                this->m_stop_worker = other.m_stop_worker;
            }
            return *this;
        }

        /**
         * @brief Allow move
         */
        CudaDeviceWorker(CudaDeviceWorker&&) = default;
        CudaDeviceWorker& operator=(CudaDeviceWorker&&) = default;


        /**
         * Starts the worker thread.
         */
        void start();

        /**
         * Stops the worker thread gracefully.
         * Ensures all pending tasks are completed before stopping.
         */
        void stop();

        /**
         * Submits a task to the worker's task queue.
         * @param task A callable object (e.g., lambda, function) to be executed by the worker.
         */
        void submit_task(std::function<void()> task);

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

        /**
         * @brief Get the ID of the CUDA device managed by this worker.
         * @return The device ID.
         */
        int get_device_id() const;

        /**
         * @brief Get the CUDA stream associated with this device worker.
         * @return The CUDA stream.
         */
        cudaStream_t get_stream() const;

        /**
         * @brief Synchronize the device.
         *
         * Waits for all operations on the device to complete.
         */
        void synchronize() const;

    private:
        /**
         * The main loop executed by the worker thread.
         * Processes tasks from the task queue until the worker is stopped.
         */
        void worker_loop();

        std::thread m_worker_thread;                      ///< The thread that runs the worker loop.
        std::mutex m_task_mutex;                          ///< Mutex for synchronizing access to the task queue.
        std::condition_variable m_task_cv;                ///< Condition variable for task queue signaling.
        std::queue<std::function<void()>> m_task_queue;   ///< FIFO queue of tasks to be executed.
        cudaStream_t m_stream;                            ///< The CUDA stream used by this worker.
        int m_device_id;                                  ///< The CUDA device ID this worker is associated with.
        bool m_stop_worker = false;                       ///< Flag to signal the worker thread to stop.
    };
}
