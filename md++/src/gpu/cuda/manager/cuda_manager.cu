
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

namespace {
    // Shared by configuration_view()'s two lookup paths: decide which
    // upload granularity covers `missing` (fields requested but not
    // already marked fresh), run it, and update the freshness bitmask.
    // FORCE/BOX have no dedicated per-field upload routine, so any
    // request touching them falls back to the coarse full copy --
    // matches today's only two existing granularities
    // (copy_to_device()/copy_pos_vel_to_device()), just chosen by
    // tracked freshness instead of a boolean picked at the call site.
    void resync_missing_fields(gpu::Configuration & mirror,
                                configuration::Configuration & conf,
                                unsigned missing) {
        if (missing == 0) return;
        if (missing & (gpu::MIRROR_FORCE | gpu::MIRROR_BOX)) {
            // Full copy_to_device() overwrites POS/VEL too -- if either
            // is currently GPU-dirty (a kernel wrote it, CPU hasn't
            // seen it yet), publish it first so this doesn't silently
            // discard that value in favour of the stale CPU copy.
            if (mirror.gpu_dirty_fields & (gpu::MIRROR_POS | gpu::MIRROR_VEL)) {
                mirror.copy_pos_vel_from_device(conf);
                mirror.gpu_dirty_fields &= ~(gpu::MIRROR_POS | gpu::MIRROR_VEL);
            }
            mirror.copy_to_device(conf);
            mirror.gpu_fresh_fields = gpu::MIRROR_ALL;
        } else {
            mirror.copy_pos_vel_to_device(conf);
            mirror.gpu_fresh_fields |= gpu::MIRROR_POS | gpu::MIRROR_VEL;
        }
    }
}

gpu::Configuration::View gpu::CudaManager::configuration_view(configuration::Configuration & conf,
                                                                unsigned read_fields) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        resync_missing_fields(*m_last_conf_gpu, conf, read_fields & ~m_last_conf_gpu->gpu_fresh_fields);
        return m_last_conf_gpu->view();
    }

    auto it = m_configurations.find(id);
    if (it == m_configurations.end()) {
        it = m_configurations.emplace(id, std::make_unique<gpu::Configuration>()).first;
        it->second->copy_to_device(conf); // full sync on first creation
        it->second->gpu_fresh_fields = gpu::MIRROR_ALL;
    } else {
        resync_missing_fields(*it->second, conf, read_fields & ~it->second->gpu_fresh_fields);
    }

    m_last_conf_id  = id;
    m_last_conf_gpu = it->second.get();
    return it->second->view();
}

void gpu::CudaManager::mark_gpu_dirty(configuration::Configuration & conf, unsigned fields) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        m_last_conf_gpu->gpu_fresh_fields |= fields;
        m_last_conf_gpu->gpu_dirty_fields |= fields;
        return;
    }

    auto it = m_configurations.find(id);
    if (it != m_configurations.end()) {
        it->second->gpu_fresh_fields |= fields;
        it->second->gpu_dirty_fields |= fields;
    }
}

void gpu::CudaManager::flush_gpu_dirty(configuration::Configuration & conf, unsigned fields) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    const unsigned to_flush = fields & mirror->gpu_dirty_fields;
    if (to_flush & (gpu::MIRROR_POS | gpu::MIRROR_VEL)) {
        mirror->copy_pos_vel_from_device(conf);
        mirror->gpu_dirty_fields &= ~(gpu::MIRROR_POS | gpu::MIRROR_VEL);
    }
    // FORCE/BOX: nothing ever marks these dirty today (mark_gpu_dirty()
    // is only called for VEL), so there's no flush routine needed for
    // them yet -- add one here if a future writer starts leaving them
    // GPU-only too.
}

void gpu::CudaManager::invalidate_gpu_mirror(configuration::Configuration & conf, unsigned fields) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        m_last_conf_gpu->gpu_fresh_fields &= ~fields;
        return;
    }

    auto it = m_configurations.find(id);
    if (it != m_configurations.end()) {
        it->second->gpu_fresh_fields &= ~fields;
    }
}

void gpu::CudaManager::sync_configuration_from_device(configuration::Configuration & conf) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        m_last_conf_gpu->copy_pos_vel_from_device(conf);
        m_last_conf_gpu->gpu_dirty_fields &= ~(gpu::MIRROR_POS | gpu::MIRROR_VEL);
        return;
    }

    auto it = m_configurations.find(id);
    if (it != m_configurations.end()) {
        it->second->copy_pos_vel_from_device(conf);
        it->second->gpu_dirty_fields &= ~(gpu::MIRROR_POS | gpu::MIRROR_VEL);
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