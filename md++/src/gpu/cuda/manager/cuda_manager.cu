
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
    // Destroy and clear every producer event guarding the bits in
    // `fields` -- called at every "this field starts fresh again" point
    // (CPU->GPU resync, zero_mirror_force(), invalidate_gpu_mirror())
    // so the per-field event lists never accumulate across steps or
    // outlive the write they were guarding.
    void clear_producer_events(gpu::Configuration & mirror, unsigned fields) {
        for (unsigned bit = 0; bit < 8; ++bit) {
            if (!(fields & (1u << bit))) continue;
            for (cudaEvent_t e : mirror.field_producer_events[bit]) cudaEventDestroy(e);
            mirror.field_producer_events[bit].clear();
        }
    }

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
        if (missing & (gpu::MIRROR_FORCE | gpu::MIRROR_BOX |
                        gpu::MIRROR_CONSTRAINT_FORCE | gpu::MIRROR_VIRIAL)) {
            // Full copy_to_device() overwrites POS/VEL too -- if either
            // is currently GPU-dirty (a kernel wrote it, CPU hasn't
            // seen it yet), publish it first so this doesn't silently
            // discard that value in favour of the stale CPU copy.
            if (mirror.gpu_dirty_fields & (gpu::MIRROR_POS | gpu::MIRROR_VEL)) {
                mirror.copy_pos_vel_from_device(conf);
                mirror.gpu_dirty_fields &= ~(gpu::MIRROR_POS | gpu::MIRROR_VEL);
            }
            mirror.copy_to_device(conf);
            // copy_to_device() also copies lattice_shifts (configuration_
            // struct.cu) even though that bit isn't part of MIRROR_ALL --
            // mark it fresh too, or the very next MIRROR_LATTICE_SHIFT
            // request would think it's missing and redundantly re-copy.
            mirror.gpu_fresh_fields = gpu::MIRROR_ALL | gpu::MIRROR_LATTICE_SHIFT;
            // Every field just got overwritten from the CPU side -- any
            // GPU producer event still pending for them is now
            // irrelevant (this resync already waited out/ superseded
            // whatever it was guarding via the cudaDeviceSynchronize()
            // inside copy_pos_vel_from_device() above, when that ran).
            clear_producer_events(mirror, gpu::MIRROR_ALL);
        } else {
            mirror.copy_pos_vel_to_device(conf);
            mirror.gpu_fresh_fields |= gpu::MIRROR_POS | gpu::MIRROR_VEL;
            clear_producer_events(mirror, gpu::MIRROR_POS | gpu::MIRROR_VEL);
        }
        // LATTICE_SHIFT is independent of the POS/VEL/coarse branches
        // above (a single persistent array, not part of either copy
        // routine) -- handled separately so requesting it never
        // triggers an unrelated full copy_to_device().
        if (missing & gpu::MIRROR_LATTICE_SHIFT) {
            mirror.copy_lattice_shifts_to_device(conf);
            mirror.gpu_fresh_fields |= gpu::MIRROR_LATTICE_SHIFT;
            clear_producer_events(mirror, gpu::MIRROR_LATTICE_SHIFT);
        }
    }

    // Insert cudaStreamWaitEvent() on `stream` for every producer event
    // guarding a field in `fields` that's already fresh (not part of
    // `missing` -- those were just handled, synchronously, by
    // resync_missing_fields() above). Pure GPU-side ordering, no CPU
    // blocking; a no-op if `stream` is null (legacy/CPU-side caller) or
    // a field has no registered producer (e.g. it became fresh via a
    // CPU upload, not a GPU kernel).
    void wait_on_fresh_producers(gpu::Configuration & mirror, unsigned read_fields,
                                  unsigned missing, cudaStream_t stream) {
        if (stream == 0) return;
        const unsigned already_fresh = read_fields & ~missing;
        for (unsigned bit = 0; bit < 8; ++bit) {
            if (!(already_fresh & (1u << bit))) continue;
            for (cudaEvent_t e : mirror.field_producer_events[bit])
                cudaStreamWaitEvent(stream, e, 0);
        }
    }
}

gpu::Configuration::View gpu::CudaManager::configuration_view(configuration::Configuration & conf,
                                                                unsigned read_fields,
                                                                cudaStream_t stream) {
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        const unsigned missing = read_fields & ~m_last_conf_gpu->gpu_fresh_fields;
        resync_missing_fields(*m_last_conf_gpu, conf, missing);
        wait_on_fresh_producers(*m_last_conf_gpu, read_fields, missing, stream);
        return m_last_conf_gpu->view();
    }

    auto it = m_configurations.find(id);
    if (it == m_configurations.end()) {
        it = m_configurations.emplace(id, std::make_unique<gpu::Configuration>()).first;
        it->second->copy_to_device(conf); // full sync on first creation
        it->second->gpu_fresh_fields = gpu::MIRROR_ALL | gpu::MIRROR_LATTICE_SHIFT;
    } else {
        const unsigned missing = read_fields & ~it->second->gpu_fresh_fields;
        resync_missing_fields(*it->second, conf, missing);
        wait_on_fresh_producers(*it->second, read_fields, missing, stream);
    }

    m_last_conf_id  = id;
    m_last_conf_gpu = it->second.get();
    return it->second->view();
}

void gpu::CudaManager::mark_gpu_dirty(configuration::Configuration & conf, unsigned fields,
                                       cudaStream_t stream) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    mirror->gpu_fresh_fields |= fields;
    mirror->gpu_dirty_fields |= fields;
    if (stream == 0) return;

    // One event per written bit (never one event shared across several
    // bits' vectors -- each vector owns and destroys its own handles,
    // so sharing would double-destroy). All recorded at the same stream
    // position, so they're functionally simultaneous.
    for (unsigned bit = 0; bit < 8; ++bit) {
        if (!(fields & (1u << bit))) continue;
        cudaEvent_t ev;
        cudaEventCreateWithFlags(&ev, cudaEventDisableTiming);
        cudaEventRecord(ev, stream);
        mirror->field_producer_events[bit].push_back(ev);
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
    if (to_flush & (gpu::MIRROR_CONSTRAINT_FORCE | gpu::MIRROR_VIRIAL)) {
        // Publishes both fields together (copy_constraint_data_from_
        // device() has no finer granularity) -- the constraint
        // algorithms always mark_gpu_dirty() both at once anyway, so
        // this never does unnecessary work in practice.
        mirror->copy_constraint_data_from_device(conf);
        mirror->gpu_dirty_fields &= ~(gpu::MIRROR_CONSTRAINT_FORCE | gpu::MIRROR_VIRIAL);
    }
    if (to_flush & gpu::MIRROR_FORCE) {
        // NonBonded/bonded-term Interactions write force directly into
        // the mirror and mark_gpu_dirty(MIRROR_FORCE) instead of each
        // keeping a private buffer and syncing it back individually --
        // publish once here, only when a downstream CPU-side algorithm
        // (or a test reading conf.current().force directly) actually
        // needs it.
        mirror->copy_forces_from_device(conf);
        mirror->gpu_dirty_fields &= ~gpu::MIRROR_FORCE;
    }
    if (to_flush & gpu::MIRROR_LATTICE_SHIFT) {
        mirror->copy_lattice_shifts_from_device(conf);
        mirror->gpu_dirty_fields &= ~gpu::MIRROR_LATTICE_SHIFT;
    }
    if (to_flush & gpu::MIRROR_ENERGY) {
        // Every GPU-native bonded/special term already atomicAdd-
        // merged its own contribution into the mirror's energy_*
        // buffers (energy_accumulate_kernels.h) -- this is the one
        // real publish point, at Energy_Calculation's own MIRROR_
        // ENERGY touch.
        mirror->copy_energy_from_device(conf);
        mirror->gpu_dirty_fields &= ~gpu::MIRROR_ENERGY;
    }
    // BOX: nothing ever marks this dirty today, so there's no flush
    // routine needed for it yet -- add one here if a future writer
    // starts leaving it GPU-only too.
}

void gpu::CudaManager::publish_cpu_virial(configuration::Configuration & conf) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    // Synchronous H2D copy (see copy_virial_to_device()'s doc comment):
    // the next GPU writer of virial_tensor (a constraint algorithm's
    // atomicAdd, issued afterward in host program order, i.e. later in
    // Algorithm_Sequence::run()) is guaranteed to be enqueued only after
    // this blocking call has already completed on-device.
    mirror->copy_virial_to_device(conf);
    // The mirror's virial_tensor now matches the CPU correction exactly
    // -- nothing to publish until a future writer marks it dirty again.
    mirror->gpu_dirty_fields &= ~gpu::MIRROR_VIRIAL;
}

void gpu::CudaManager::invalidate_gpu_mirror(configuration::Configuration & conf, unsigned fields) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    mirror->gpu_fresh_fields &= ~fields;
    // No longer fresh -- any producer event still pending for these
    // bits is stale; the next writer re-populates when it calls
    // mark_gpu_dirty() again.
    clear_producer_events(*mirror, fields);
}

void gpu::CudaManager::clear_stale_producer_events(configuration::Configuration & conf, unsigned fields) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    clear_producer_events(*mirror, fields);
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

void gpu::CudaManager::exchange_mirror_state(configuration::Configuration & conf) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return; // no mirror yet -- nothing to swap

    mirror->exchange_state();
    // The producer-event lists are indexed by field bit, not by
    // current/old -- current.force's events described current's
    // contents, which is now old's. Field_producer_events tracks
    // "who last wrote MIRROR_FORCE" regardless of which half that was,
    // so no bookkeeping change is needed here beyond the struct swap
    // itself: a consumer requesting MIRROR_FORCE still finds the right
    // (now current) events pending, since gpu_fresh_fields/
    // gpu_dirty_fields and field_producer_events are per-mirror, not
    // per-half.
}

void gpu::CudaManager::zero_mirror_force(configuration::Configuration & conf) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return; // no mirror yet -- first configuration_view() call will build one from the already-zeroed CPU force

    if (mirror->current.force.size() > 0) {
        cudaMemset(mirror->current.force.data(), 0,
                   mirror->current.force.size() * sizeof(FPH3_TYPE));
    }
    if (mirror->current.virial_tensor) {
        cudaMemset(mirror->current.virial_tensor, 0, sizeof(FPH9_TYPE));
    }
    // Zeroed here (not per-constraint-algorithm) for the same reason as
    // force/virial above: called once per step, before the swap that
    // makes this "current" half into "old" -- exactly the half the
    // constraint algorithms (which run after Leap_Frog_Velocity's swap)
    // publish their own constraint_force contribution into. Solute vs
    // solvent constraint algorithms write disjoint atom ranges (exactly
    // one of SHAKE/LINCS is active for solute, one of M_SHAKE/SETTLE for
    // solvent), so no cross-algorithm accumulation ordering is needed --
    // each just needs its own range zeroed first.
    if (mirror->current.constraint_force.size() > 0) {
        cudaMemset(mirror->current.constraint_force.data(), 0,
                   mirror->current.constraint_force.size() * sizeof(FPH3_TYPE));
    }
    // Fresh step: every producer event guarding last step's FORCE value
    // is now meaningless (the buffer was just zeroed by this call, on
    // the default stream -- see this method's ordering note below).
    // Each Interaction's own mark_gpu_dirty(MIRROR_FORCE, its_stream)
    // this step re-populates it.
    clear_producer_events(*mirror, gpu::MIRROR_FORCE);
    clear_producer_events(*mirror, gpu::MIRROR_CONSTRAINT_FORCE);
    // Same reasoning applies to MIRROR_VIRIAL: every bonded term/
    // NonBonded/active constraint algorithm calls mark_gpu_dirty(...,
    // MIRROR_VIRIAL, its_stream) once per step (virial_accumulate_
    // kernels.h), each appending one more event -- without this clear,
    // field_producer_events[VIRIAL] grows by ~6 events every step
    // forever (real cudaEventCreateWithFlags()/cudaEventRecord()
    // resources, never destroyed), the same O(steps^2)-event-leak
    // pattern MIRROR_LATTICE_SHIFT was already fixed for elsewhere in
    // this file -- missed here when virial_tensor's on-device
    // accumulation was added. (The dominant cost turned out to be a
    // sibling leak on MIRROR_POS/MIRROR_VEL, fixed in lattice_shift_
    // gpu.cc/leap_frog_gpu.cc -- this one alone is a smaller, real fix,
    // not the full story.)
    clear_producer_events(*mirror, gpu::MIRROR_VIRIAL);
}

void gpu::CudaManager::ensure_energy_groups(configuration::Configuration & conf, unsigned num_groups) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return;

    mirror->resize_energy_groups(num_groups);
}

gpu::EnergyMirrorPtrs gpu::CudaManager::energy_mirror_ptrs(configuration::Configuration & conf) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror) return gpu::EnergyMirrorPtrs{};

    return gpu::EnergyMirrorPtrs{
        mirror->energy_bond, mirror->energy_angle, mirror->energy_improper,
        mirror->energy_dihedral, mirror->energy_posrest};
}

void gpu::CudaManager::zero_mirror_energy(configuration::Configuration & conf) {
    gpu::Configuration * mirror = nullptr;
    const std::size_t id = conf.id();

    if (id == m_last_conf_id && m_last_conf_gpu) {
        mirror = m_last_conf_gpu;
    } else {
        auto it = m_configurations.find(id);
        if (it != m_configurations.end()) mirror = it->second.get();
    }
    if (!mirror || mirror->energy_num_groups == 0) return;

    mirror->zero_energy(0);
    // Fresh step: every producer event guarding last step's energy_*
    // values is now meaningless -- same reasoning as zero_mirror_
    // force()'s MIRROR_VIRIAL clear.
    clear_producer_events(*mirror, gpu::MIRROR_ENERGY);
}

int * gpu::CudaManager::constraint_error_flag_slot(unsigned slot) {
    if (m_constraint_error_flags.size() == 0) {
        m_constraint_error_flags.resize(gpu::NUM_CONSTRAINT_ERROR_SLOTS);
        cudaMemset(m_constraint_error_flags.data(), 0,
                   gpu::NUM_CONSTRAINT_ERROR_SLOTS * sizeof(int));
    }
    return m_constraint_error_flags.data() + slot;
}

void gpu::CudaManager::zero_constraint_error_flags() {
    if (m_constraint_error_flags.size() == 0) return; // nothing allocated yet
    cudaMemset(m_constraint_error_flags.data(), 0,
               gpu::NUM_CONSTRAINT_ERROR_SLOTS * sizeof(int));
}

bool gpu::CudaManager::check_constraint_error_flags(std::vector<int> & out_codes) {
    if (m_constraint_error_flags.size() == 0) return false; // no GPU constraint ran this build/run

    // The one sync for every constraint algorithm's deferred status,
    // instead of each checking (and syncing on) its own flag right
    // after its own kernels.
    cudaDeviceSynchronize();

    out_codes.assign(gpu::NUM_CONSTRAINT_ERROR_SLOTS, 0);
    for (unsigned s = 0; s < gpu::NUM_CONSTRAINT_ERROR_SLOTS; ++s) {
        out_codes[s] = m_constraint_error_flags[s];
    }

    bool fatal = false;
    for (unsigned s = 0; s < gpu::NUM_CONSTRAINT_ERROR_SLOTS; ++s) {
        if (out_codes[s] != 0 && gpu::constraint_error_slot_is_fatal(s)) fatal = true;
    }
    return fatal;
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