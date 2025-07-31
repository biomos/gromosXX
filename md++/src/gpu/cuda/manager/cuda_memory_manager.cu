
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "gpu.h"

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"
#include "cuda_memory_manager.h"

gpu::CudaMemoryManager::CudaMemoryManager(int device_id) {}

gpu::CudaMemoryManager::~CudaMemoryManager() {
    cudaSetDevice(m_device_id);
    // Free all tracked allocations to avoid leaks if not freed yet
    std::lock_guard<std::mutex> lock(m_mutex);
    for (const auto&  alloc_pair : m_allocations) {
      cudaFree(alloc_pair.first);
    }
    m_allocations.clear();
}

template <typename T>
gpu::CudaMemoryManager::cumm_unique_ptr<T> gpu::CudaMemoryManager::allocate(size_t count) {
    cudaSetDevice(m_device_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    T* devptr;
    cudaMalloc(&devptr, sizeof(T));
    m_allocations[devptr] = sizeof(T);
    return m_unique_ptr(devptr);
};


template<typename T>
gpu::CudaMemoryManager::cumm_unique_ptr<std::remove_cv_t<T>> gpu::CudaMemoryManager::copy_to_device(const T& src) {
    using rawT = std::remove_cv_t<T>;
    cumm_unique_ptr<rawT> device_ptr = allocate<rawT>();
    std::lock_guard<std::mutex> lock(m_mutex);
    cudaSetDevice(m_device_id);
    cudaMemcpy(device_ptr.get(), &src, sizeof(rawT), cudaMemcpyHostToDevice);
    return device_ptr;
}

template<typename C, typename T, typename U>
gpu::CudaMemoryManager::cumm_unique_ptr<U> gpu::CudaMemoryManager::copy_to_device(const C& vec) {
    cudaSetDevice(m_device_id);
    cumm_unique_ptr<U> device_ptr = allocate<U>(vec.size());
    std::lock_guard<std::mutex> lock(m_mutex);

    if constexpr (!std::is_same<T, U>::value) {
        std::vector<U> tmp(vec.size());
        std::transform(vec.begin(), vec.end(), tmp.begin(),
                    [](const T& t) { return static_cast<U>(t); });

        cudaMemcpy(device_ptr.get(), tmp.data(), sizeof(U) * tmp.size(), cudaMemcpyHostToDevice);
    } else {
        cudaMemcpy(device_ptr.get(), vec.data(), sizeof(U) * vec.size(), cudaMemcpyHostToDevice);
    }
    return device_ptr;
}

template <typename T>
T* gpu::CudaMemoryManager::allocate_managed(const T& host_data) {
    cudaSetDevice(m_device_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    T* devptr;
    cudaMallocManaged(&devptr, sizeof(T));
    *devptr = host_data;
    m_allocations[devptr] = sizeof(T);
    return devptr;
};

void gpu::CudaMemoryManager::deallocate(void* device_ptr) {
    std::lock_guard<std::mutex> lock(m_mutex);

    auto it = m_allocations.find(device_ptr);
    if (it == m_allocations.end()) {
        throw std::runtime_error("Attempt to free memory that was not allocated by this CudaMemoryManager.");
    }

    cudaFree(device_ptr);
    m_allocations.erase(it);
}

void gpu::CudaMemoryManager::validate(void* ptr) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_allocations.find(ptr) == m_allocations.end())
        throw std::runtime_error("Pointer not owned by this CudaMemoryManager");
}

void gpu::CudaMemoryManager::print_allocations() const {
    std::lock_guard<std::mutex> lock(m_mutex);

    for (const auto& [ptr, size] : m_allocations) {
        std::cout << "Device pointer: " << ptr << ", Size: " << size << " bytes\n";
    }
}

// void gpu::CudaMemoryManager::init() {}

// void* gpu::CudaMemoryManager::allocate_device_memory(size_t size) {
//     void* device_ptr = nullptr;
//     cudaSetDevice(m_device_id);
//     CUDA_CHECK(cudaMalloc(&device_ptr, size));
//     return device_ptr;
// }

// void gpu::CudaMemoryManager::free_device_memory(void* device_ptr) {
//     CUDA_CHECK(cudaFree(device_ptr));
// }

void* gpu::CudaMemoryManager::allocate_pinned_memory(size_t size) {
    void* host_ptr = nullptr;
    CUDA_CHECK(cudaHostAlloc(&host_ptr, size, cudaHostAllocDefault));
    return host_ptr;
}

void gpu::CudaMemoryManager::free_pinned_memory(void* host_ptr) {
    CUDA_CHECK(cudaFreeHost(host_ptr));
}

void gpu::CudaMemoryManager::copy_to_device(void* device_ptr, const void* host_ptr, size_t size, CUSTREAM stream) {
    CUDA_CHECK(cudaMemcpy(device_ptr, host_ptr, size, cudaMemcpyHostToDevice));
}

void gpu::CudaMemoryManager::copy_to_host(void* host_ptr, const void* device_ptr, size_t size, CUSTREAM stream) {
    CUDA_CHECK(cudaMemcpy(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost));
}

void gpu::CudaMemoryManager::copy_device_to_device(void* dest_device_ptr, const void* src_device_ptr, size_t size) {
    CUDA_CHECK(cudaMemcpy(dest_device_ptr, src_device_ptr, size, cudaMemcpyDeviceToDevice));
}

void gpu::CudaMemoryManager::async_copy_to_device(void* device_ptr, const void* host_ptr, size_t size, cudaStream_t stream) {
    CUDA_CHECK(cudaMemcpyAsync(device_ptr, host_ptr, size, cudaMemcpyHostToDevice, stream));
}

void gpu::CudaMemoryManager::async_copy_to_host(void* host_ptr, const void* device_ptr, size_t size, cudaStream_t stream) {
    CUDA_CHECK(cudaMemcpyAsync(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost, stream));
}

void gpu::CudaMemoryManager::query_memory(size_t& free_memory, size_t& total_memory) const {
    size_t free_mem = 0;
    size_t total_mem = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    free_memory = free_mem;
    total_memory = total_mem;
}
