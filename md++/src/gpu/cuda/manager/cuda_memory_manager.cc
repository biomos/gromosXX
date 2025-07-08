#include "cuda_memory_manager.h"
#include <stdexcept>
#include <sstream>

gpu::CudaMemoryManager::CudaMemoryManager() {
    // Constructor: No specific initialization required for now
}

gpu::CudaMemoryManager::~CudaMemoryManager() {
    // Destructor: No specific cleanup required for now
}

void* gpu::CudaMemoryManager::allocate_device_memory(size_t size) {
    void* device_ptr = nullptr;
    cudaError_t result = cudaMalloc(&device_ptr, size);
    check_cuda_error(result, "Failed to allocate device memory.");
    return device_ptr;
}

void gpu::CudaMemoryManager::free_device_memory(void* device_ptr) {
    cudaError_t result = cudaFree(device_ptr);
    check_cuda_error(result, "Failed to free device memory.");
}

void* gpu::CudaMemoryManager::allocate_pinned_memory(size_t size) {
    void* host_ptr = nullptr;
    cudaError_t result = cudaHostAlloc(&host_ptr, size, cudaHostAllocDefault);
    check_cuda_error(result, "Failed to allocate pinned host memory.");
    return host_ptr;
}

void gpu::CudaMemoryManager::free_pinned_memory(void* host_ptr) {
    cudaError_t result = cudaFreeHost(host_ptr);
    check_cuda_error(result, "Failed to free pinned host memory.");
}

void gpu::CudaMemoryManager::copy_to_device(void* device_ptr, const void* host_ptr, size_t size) {
    cudaError_t result = cudaMemcpy(device_ptr, host_ptr, size, cudaMemcpyHostToDevice);
    check_cuda_error(result, "Failed to copy data from host to device.");
}

void gpu::CudaMemoryManager::copy_to_host(void* host_ptr, const void* device_ptr, size_t size) {
    cudaError_t result = cudaMemcpy(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost);
    check_cuda_error(result, "Failed to copy data from device to host.");
}

void gpu::CudaMemoryManager::copy_device_to_device(void* dest_device_ptr, const void* src_device_ptr, size_t size) {
    cudaError_t result = cudaMemcpy(dest_device_ptr, src_device_ptr, size, cudaMemcpyDeviceToDevice);
    check_cuda_error(result, "Failed to copy data between device memory locations.");
}

void gpu::CudaMemoryManager::async_copy_to_device(void* device_ptr, const void* host_ptr, size_t size, cudaStream_t stream) {
    cudaError_t result = cudaMemcpyAsync(device_ptr, host_ptr, size, cudaMemcpyHostToDevice, stream);
    check_cuda_error(result, "Failed to perform asynchronous copy from host to device.");
}

void gpu::CudaMemoryManager::async_copy_to_host(void* host_ptr, const void* device_ptr, size_t size, cudaStream_t stream) {
    cudaError_t result = cudaMemcpyAsync(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost, stream);
    check_cuda_error(result, "Failed to perform asynchronous copy from device to host.");
}

void gpu::CudaMemoryManager::query_memory(size_t& free_memory, size_t& total_memory) const {
    size_t free_mem = 0;
    size_t total_mem = 0;
    cudaError_t result = cudaMemGetInfo(&free_mem, &total_mem);
    check_cuda_error(result, "Failed to query device memory.");
    free_memory = free_mem;
    total_memory = total_mem;
}

void gpu::CudaMemoryManager::check_cuda_error(cudaError_t result, const std::string& message) const {
    if (result != cudaSuccess) {
        std::ostringstream oss;
        oss << message << " CUDA Error: " << cudaGetErrorString(result);
        throw std::runtime_error(oss.str());
    }
}