
#include <stdexcept>
#include <sstream>

#include "gpu.h"

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"
#include "cuda_memory_manager.h"

gpu::CudaMemoryManager::CudaMemoryManager(int device_id) {}

gpu::CudaMemoryManager::~CudaMemoryManager() {}

void gpu::CudaMemoryManager::init() {
    pos_umem.resize(10000,0);
    pos_host.resize(10000,0);
}

void* gpu::CudaMemoryManager::allocate_device_memory(size_t size) {
    void* device_ptr = nullptr;
    cudaSetDevice(m_device_id);
    CUDA_CHECK(cudaMalloc(&device_ptr, size));
    return device_ptr;
}

void gpu::CudaMemoryManager::free_device_memory(void* device_ptr) {
    CUDA_CHECK(cudaFree(device_ptr));
}

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
