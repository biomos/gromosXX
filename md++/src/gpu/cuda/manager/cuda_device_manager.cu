
#include <sstream>

#include "gpu.h"

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"
#include "cuda_device_manager.h"

gpu::CudaDeviceManager::CudaDeviceManager(int device_id) : m_device_id(device_id) {
    // setting device id is not neccessary as we do not interact with the device yet
    // CUDA_CHECK(cudaSetDevice(m_device_id));
    // add one worker
    add_worker();
}

gpu::CudaDeviceManager::~CudaDeviceManager() {
    // No explicit cleanup needed for device manager
    //cudaDeviceReset();
}

cudaDeviceProp gpu::CudaDeviceManager::get_device_properties() const {
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, m_device_id);
    return properties;
}

/**
 * Adds new worker to the vector of workers
 */
void gpu::CudaDeviceManager::add_worker() {
    std::lock_guard<std::mutex> lock(m_workers_mutex);
    m_workers.emplace_back(CudaDeviceWorker(m_device_id));
}

/**
 * Remove worker from the vector of workers
 */
void gpu::CudaDeviceManager::remove_worker() {
    std::lock_guard<std::mutex> lock(m_workers_mutex);
    if (m_workers.empty()) {
        throw std::runtime_error("No workers to remove.");
    }

    // Stop the last worker
    m_workers.back().stop();

    // Remove the last worker from the vector
    m_workers.pop_back();
}

/**
 * Submits a task to the next worker
 */
void gpu::CudaDeviceManager::submit_task(std::function<void()> task) {
    size_t worker_index;
    {
        std::lock_guard<std::mutex> lock(m_workers_mutex); // Protect access to next_worker
        worker_index = m_next_worker;
        m_next_worker = (m_next_worker + 1) % m_workers.size();
    }
    m_workers[worker_index].submit_task(std::move(task));
}

int gpu::CudaDeviceManager::get_device_id() const {
    return m_device_id;
}

// cudaStream_t gpu::CudaDeviceManager::create_stream() {
//     cudaStream_t stream;
//     cudaStreamCreate(&stream);
//     return stream;
// }

// void gpu::CudaDeviceManager::destroy_stream(cudaStream_t stream) {
//     cudaStreamDestroy(stream);
// }

std::string gpu::CudaDeviceManager::get_device_description() const {
    cudaDeviceProp properties = get_device_properties();

    std::ostringstream oss;
    oss << "Device " << m_device_id << ": " << properties.name << "\n"
        << "  Compute Capability: " << properties.major << "." << properties.minor << "\n"
        << "  Total Global Memory: " << (properties.totalGlobalMem / (1024 * 1024)) << " MB\n"
        << "  Multiprocessors: " << properties.multiProcessorCount << "\n"
        << "  Max Threads per Block: " << properties.maxThreadsPerBlock << "\n"
        << "  Max Threads per SM: " << properties.maxThreadsPerMultiProcessor << "\n";
    return oss.str();
}