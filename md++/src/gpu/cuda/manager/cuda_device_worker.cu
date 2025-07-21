
#include "gpu.h"

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"
#include "cuda_device_worker.h"

gpu::CudaDeviceWorker::CudaDeviceWorker(int device_id) :
                        m_stream(nullptr), m_device_id(device_id), m_stop_worker(false) {
    // Set the active device
    CUDA_CHECK(cudaSetDevice(m_device_id));

    // Create a CUDA stream
    CUDA_CHECK(cudaStreamCreate(&m_stream));
}

gpu::CudaDeviceWorker::~CudaDeviceWorker() {
    // Destroy the CUDA stream
    if (m_stream) {
        CUDA_CHECK(cudaStreamDestroy(m_stream));
    }
}

/**
 * Starts the worker thread.
 */
void gpu::CudaDeviceWorker::start() {
    m_worker_thread = std::thread(&CudaDeviceWorker::worker_loop, this);
}

/**
 * Stops the worker thread gracefully.
 * Ensures all pending tasks are completed before stopping.
 */
void gpu::CudaDeviceWorker::stop() {
    {
        std::lock_guard<std::mutex> lock(m_task_mutex);
        m_stop_worker = true;
    }
    m_task_cv.notify_all();
    if (m_worker_thread.joinable()) {
        m_worker_thread.join();
    }
}

void gpu::CudaDeviceWorker::submit_task(std::function<void()> task) {
    {
        std::lock_guard<std::mutex> lock(m_task_mutex);
        m_task_queue.push(std::move(task));
    }
    m_task_cv.notify_one();
}

/**
 * Worker loop for processing tasks
 */
void gpu::CudaDeviceWorker::worker_loop() {
    CUDA_CHECK(cudaSetDevice(m_device_id));
    while (true) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock(m_task_mutex);

            // Wait until there is a task or m_stop_worker is true
            m_task_cv.wait(lock, [this]() { return !m_task_queue.empty() || m_stop_worker; });

            // If m_stop_worker is true and no tasks remain, exit the loop
            if (m_stop_worker && m_task_queue.empty()) {
                break;
            }

            // Get the next task
            if (!m_task_queue.empty()) {
                task = std::move(m_task_queue.front());
                m_task_queue.pop();
            }
        }

        // Execute the task outside the lock
        if (task) {
            task();
        }
    }
}

int gpu::CudaDeviceWorker::get_device_id() const {
    return m_device_id;
}

cudaStream_t gpu::CudaDeviceWorker::get_stream() const {
    return m_stream;
}

void gpu::CudaDeviceWorker::synchronize() const {
    CUDA_CHECK(cudaDeviceSynchronize());
}

template <typename KernelFunc, typename... Args>
void gpu::CudaDeviceWorker::launch_kernel(KernelFunc kernel, dim3 grid_dim, dim3 block_dim, Args... args, size_t shared_mem_size) {
    kernel<<<grid_dim, block_dim, shared_mem_size, m_stream>>>(args...);

    // Check for errors after kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(err));
    }
}
