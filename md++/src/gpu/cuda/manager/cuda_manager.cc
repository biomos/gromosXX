
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

// Genuine no-ops, not DISABLED_VOID() -- these two are called
// unconditionally, on every algorithm, every step, by
// Algorithm_Sequence::run() regardless of accelerator or whether this
// build even has CUDA compiled in (see the doc comments in
// cuda_manager.h). Matches the CUDA-enabled implementation's own
// runtime behavior when a Configuration has no GPU mirror registered
// (cuda_manager.cu: "if (!mirror) return;") -- there's simply nothing
// to flush/invalidate in a CPU-only build, not an attempted use of
// unavailable CUDA functionality. Using DISABLED_VOID() here used to
// log a "CUDA is disabled" critical message on every single call --
// thousands of them over an ordinary CPU-only run -- found via a real
// end-to-end run (extended_test/).
void gpu::CudaManager::flush_gpu_dirty(configuration::Configuration & conf, unsigned fields) {
}

void gpu::CudaManager::publish_cpu_virial(configuration::Configuration & conf) {
}

void gpu::CudaManager::invalidate_gpu_mirror(configuration::Configuration & conf, unsigned fields) {
}

void gpu::CudaManager::clear_stale_producer_events(configuration::Configuration & conf, unsigned fields) {
}

void gpu::CudaManager::exchange_mirror_state(configuration::Configuration & conf) {
}

// Same "genuine no-op, called unconditionally every step" reasoning as
// flush_gpu_dirty()/invalidate_gpu_mirror() above -- Algorithm_Sequence
// ::run() calls these regardless of accelerator.
void gpu::CudaManager::zero_mirror_force(configuration::Configuration & conf) {
}

void gpu::CudaManager::zero_constraint_error_flags() {
}

bool gpu::CudaManager::check_constraint_error_flags(std::vector<int> & out_codes) {
    return false;
}