
#include <memory>
#include <stdexcept>
#include <sstream>

#include "gpu/cuda/cuheader.h"
#include "cuda_manager.h"

#ifdef USE_CUDA
bool gpu::CudaManager::is_enabled_ = false;
#endif

gpu::CudaManager::CudaManager()
    : device_manager_(std::make_unique<CudaDeviceManager>()), memory_manager_() {
#ifndef USE_CUDA
        DISABLED_VOID();
#endif
    }

gpu::CudaManager::~CudaManager() {
    // Ensure all device workers are cleaned up
    device_workers_.clear();
}

void gpu::CudaManager::init(const std::vector<int>& device_ids) {
    // Query available devices
    int available_device_count = device_manager_->get_device_count();
    if (available_device_count == 0) {
        io::messages.add("No CUDA devices available.",
            "CudaManager", io::message::error);
    }

    // Determine which devices to initialize
    std::vector<int> devices_to_initialize = device_ids.empty()
        ? std::vector<int>(available_device_count)
        : device_ids;

    if (device_ids.empty()) {
        for (int i = 0; i < available_device_count; ++i) {
            devices_to_initialize[i] = i;
        }
    }

    // Initialize workers for each device
    for (int device_id : devices_to_initialize) {
        // validate_device_id(device_id);
        device_manager_->set_device(device_id);
        device_workers_[device_id] = std::make_unique<CudaDeviceWorker>(device_id);
    }

    memory_manager_.init();
}

size_t gpu::CudaManager::get_device_count() const {
    return device_workers_.size();
}

CUSTREAM gpu::CudaManager::get_stream(int device_id) const {
    validate_device_id(device_id);
    return device_workers_.at(device_id)->get_stream();
}

void gpu::CudaManager::synchronize_all() {
    for (const auto& [device_id, worker] : device_workers_) {
        worker->synchronize();
    }
}

void gpu::CudaManager::synchronize_device(int device_id) {
    validate_device_id(device_id);
    device_workers_.at(device_id)->synchronize();
}

std::vector<std::string> gpu::CudaManager::get_active_device_descriptions() const {
    std::vector<std::string> descriptions;
    for (const auto& [device_id, worker] : device_workers_) {
        descriptions.push_back(device_manager_->get_device_description(device_id));
    }
    return descriptions;
}

void gpu::CudaManager::validate_device_id(int device_id) const {
    if (device_workers_.find(device_id) == device_workers_.end()) {
        throw std::invalid_argument("Invalid device ID: " + std::to_string(device_id));
    }
}