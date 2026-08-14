
#include "stdheader.h"

#include "algorithm/algorithm.h"

// topology.h/configuration.h no longer pull in any gpu/cuda/... headers
// (PLAN.md §3.2 removed that circular dependency), so including the full
// definitions here -- needed for topo.id()/conf.id() and to pass topo/conf
// to gpu::Topology/gpu::Configuration's own constructors/update methods --
// is safe. Included before any gpu/cuda/... header (below): configuration.h
// pulls in <complex> (via mesh.h -> math/fft.h), which nvcc mis-parses if
// something in the gpu/cuda math headers (device sin/cos/sqrt overloads)
// has already been seen first in this translation unit -- same ordering
// every other .cu file in this tree that includes both already uses.
#include "topology/topology.h"
#include "configuration/configuration.h"

#include <memory>
#include <stdexcept>
#include <sstream>

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/utils.h"

#include "cuda_device_manager.h"
#include "cuda_device_worker.h"
#include "cuda_memory_manager.h"

#include "cuda_manager.h"

#include "cuda_manager.tcc" // Include template implementations

gpu::CudaManager::CudaManager() {}

/**
 * @brief Allow shallow copy constructor, but warn
 */
gpu::CudaManager::CudaManager(const gpu::CudaManager& other) {
    // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
    std::cerr << "Warning: Shallow copy of CudaManager at " << __FILE__
            << ":" << __LINE__ << " in function " << __func__ << std::endl;
    this->m_device_managers = other.m_device_managers;
}

/**
 * @brief Allow shallow assignment operator, but warn
 */
gpu::CudaManager& gpu::CudaManager::operator=(const gpu::CudaManager& other) {
    // CUDA_MANAGER_COPY_WARNING; // Compile-time warning
    if (this != &other) {
        std::cerr << "Warning: Shallow copy assignment of CudaManager at " << __FILE__
                << ":" << __LINE__ << " in function " << __func__ << std::endl;
        // Perform shallow copy
        this->m_device_managers = other.m_device_managers;
    }
    return *this;
}

gpu::Topology::View gpu::CudaManager::topology_view(const topology::Topology & topo,
                                                     bool force_resync) {
    const std::size_t id = topo.id();

    if (id == m_last_topo_id && m_last_topo_gpu) {
        if (force_resync) m_last_topo_gpu->update(topo);
        return m_last_topo_gpu->view();
    }

    auto it = m_topologies.find(id);
    if (it == m_topologies.end()) {
        it = m_topologies.emplace(id, std::make_unique<gpu::Topology>(topo)).first;
    } else if (force_resync) {
        it->second->update(topo);
    }

    m_last_topo_id  = id;
    m_last_topo_gpu = it->second.get();
    return it->second->view();
}

gpu::Configuration::View gpu::CudaManager::configuration_view(configuration::Configuration & conf,
                                                                bool sync_pos_vel,
                                                                bool full_resync) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        if (full_resync) m_last_conf_gpu->copy_to_device(conf);
        else if (sync_pos_vel) m_last_conf_gpu->copy_pos_vel_to_device(conf);
        return m_last_conf_gpu->view();
    }

    auto it = m_configurations.find(id);
    if (it == m_configurations.end()) {
        it = m_configurations.emplace(id, std::make_unique<gpu::Configuration>()).first;
        it->second->copy_to_device(conf); // full sync on first creation
    } else if (full_resync) {
        it->second->copy_to_device(conf);
    } else if (sync_pos_vel) {
        it->second->copy_pos_vel_to_device(conf);
    }

    m_last_conf_id  = id;
    m_last_conf_gpu = it->second.get();
    return it->second->view();
}

void gpu::CudaManager::sync_configuration_from_device(configuration::Configuration & conf) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        m_last_conf_gpu->copy_pos_vel_from_device(conf);
        return;
    }

    auto it = m_configurations.find(id);
    if (it != m_configurations.end()) {
        it->second->copy_pos_vel_from_device(conf);
    }
}

void gpu::CudaManager::init(const std::vector<int>& device_ids) {
    int deviceCount = 0;
    CUDA_CHECK(cudaGetDeviceCount(&deviceCount));
    // Query available devices
    if (deviceCount == 0) {
        io::messages.add("No CUDA devices available.",
            "CudaManager", io::message::error);
    }

    // Determine which devices to initialize
    std::vector<int> devices_to_initialize = device_ids.empty()
        ? std::vector<int>(deviceCount)
        : device_ids;

    if (device_ids.empty()) {
        for (int i = 0; i < deviceCount; ++i) {
            devices_to_initialize[i] = i;
        }
    }

    // Initialize workers for each device
    for (int device_id : devices_to_initialize) {
        // validate_device_id(device_id);
        // auto [it, inserted] = m_device_managers.emplace(device_id, std::make_unique<gpu::CudaDeviceManager>(device_id));
        // if (!inserted) {
        //     throw std::runtime_error("Duplicate device ID: " + std::to_string(device_id));
        // }
    }

    // memory_manager_.init();
}

size_t gpu::CudaManager::get_device_count() const {
    return m_device_managers.size();
}

// gpu::CUSTREAM gpu::CudaManager::get_stream(int device_id) const {
//     validate_device_id(device_id);
//     return m_device_managers.at(device_id).get_stream();
// }

void gpu::CudaManager::synchronize_all() {
    for (const auto& [device_id, device_manager] : m_device_managers) {
        device_manager->synchronize();
    }
}

void gpu::CudaManager::synchronize_device(int device_id) {
    validate_device_id(device_id);
    m_device_managers.at(device_id)->synchronize();
}

std::vector<std::string> gpu::CudaManager::get_active_device_descriptions() const {
    std::vector<std::string> descriptions;
    for (const auto& [device_id, device_manager] : m_device_managers) {
        descriptions.push_back(device_manager->get_device_description());
    }
    return descriptions;
}

void gpu::CudaManager::validate_device_id(int device_id) const {
    if (m_device_managers.find(device_id) == m_device_managers.end()) {
        throw std::invalid_argument("Invalid device ID: " + std::to_string(device_id));
    }
}