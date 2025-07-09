

#include "gpu/cuda/cuheader.h"
#include "cuda_memory_manager.h"

gpu::CudaMemoryManager::CudaMemoryManager() {}

gpu::CudaMemoryManager::~CudaMemoryManager() {}

void* gpu::CudaMemoryManager::allocate_device_memory(size_t size) {
    return DISABLED(void*);
}

void gpu::CudaMemoryManager::free_device_memory(void* device_ptr) {
    DISABLED_VOID();
}

void* gpu::CudaMemoryManager::allocate_pinned_memory(size_t size) {
    return DISABLED(void*);
}

void gpu::CudaMemoryManager::free_pinned_memory(void* host_ptr) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::copy_to_device(void* device_ptr, const void* host_ptr, size_t size, CUSTREAM stream) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::copy_to_host(void* host_ptr, const void* device_ptr, size_t size, CUSTREAM stream) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::copy_device_to_device(void* dest_device_ptr, const void* src_device_ptr, size_t size) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::async_copy_to_device(void* device_ptr, const void* host_ptr, size_t size, CUSTREAM stream) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::async_copy_to_host(void* host_ptr, const void* device_ptr, size_t size, CUSTREAM stream) {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::query_memory(size_t& free_memory, size_t& total_memory) const {
    DISABLED_VOID();
}

void gpu::CudaMemoryManager::check_cuda_error(CUERROR result, const std::string& message) const {
    DISABLED_VOID();
}