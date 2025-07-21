
#include <memory>
#include <stdexcept>
#include <sstream>

#include "gpu/cuda/cuheader.h"
#include "cuda_manager.h"

gpu::CudaManager::CudaManager() {}

/**
 * @brief Allow shallow copy constructor, but warn
 */
gpu::CudaManager::CudaManager(const gpu::CudaManager& other) {}

/**
 * @brief Allow shallow assignment operator, but warn
 */
gpu::CudaManager& gpu::CudaManager::operator=(const gpu::CudaManager& other) {
    return *this;
}

// gpu::CudaManager::~CudaManager() {}

void gpu::CudaManager::init(const std::vector<int>& device_ids) {
    DISABLED_VOID();
}

size_t gpu::CudaManager::get_device_count() const {
    return DISABLED(size_t);
}

// gpu::CUSTREAM gpu::CudaManager::get_stream(int device_id) const {
//     return DISABLED(gpu::CUSTREAM);
// }

void gpu::CudaManager::synchronize_all() {
    DISABLED_VOID();
}

void gpu::CudaManager::synchronize_device(int device_id) {
    DISABLED_VOID();
}

std::vector<std::string> gpu::CudaManager::get_active_device_descriptions() const {
    return DISABLED(std::vector<std::string>);
}

void gpu::CudaManager::validate_device_id(int device_id) const {
    DISABLED_VOID();
}