#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include "gpu/cuda/cuheader.h"
#include "gpu/mirror_fields.h"
#include "gpu/constraint_error_slots.h"

#ifdef USE_CUDA
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/memory/cuvector.h"
#endif

namespace topology {
    class Topology;
}
namespace configuration {
    class Configuration;
}

namespace gpu {
    /**
     * @brief Raw device pointers to gpu::Configuration's per-energy-
     * group scratch buffers (see that class's doc comment) -- plain
     * pointers, not a View wrapper, since callers only ever pass them
     * straight into energy_accumulate_kernels.h's launch_accumulate_
     * energy(), never index them from a __global__ kernel through this
     * struct itself. All null if `conf` has no mirror yet.
     */
    struct EnergyMirrorPtrs {
        double* bond = nullptr;
        double* angle = nullptr;
        double* improper = nullptr;
        double* dihedral = nullptr;
        double* posrest = nullptr;
        // Flattened [gi * num_groups + gj], see gpu::Configuration::
        // energy_lj/energy_crf's doc comment.
        double* lj = nullptr;
        double* crf = nullptr;
    };
}

#define CUDA_VARIABLE_DISABLED() disabled(__FILE__, __LINE__, __func__)
namespace gpu {
    class CudaDeviceManager;
    // class CudaMemoryManager;
    // class CudaDeviceWorker;
    /**
     * @class CudaManager
     * @brief High-level orchestrator for multi-GPU CUDA operations in a molecular dynamics program.
     *
     * The CudaManager class initializes the CUDA environment, manages multiple devices, and coordinates
     * memory and kernel execution across GPUs. It provides a clean interface for interacting with CUDA
     * resources and abstracts away low-level details.
     * 
     * This is the only interface for host code to access CUDA, to encapsulate all cuda code separately
     * and not expose it to the CPU-only code.
     * Function calls, that submit to GPU should go only and only over this manager.
     * Probably should be called CudaInterface, we will see.
     * 
     */

    #define CUDA_MANAGER_COPY_WARNING \
    _Pragma("message(\"Warning: Shallow copy of CudaManager detected.\")")

    class CudaManager {
        public:
            /**
             * @brief Constructor
             */
            CudaManager();

            /**
             * Destructor
             */
            // ~CudaManager();

            /**
             * @brief Disable copy
             */
            // CudaManager(const CudaManager&) = delete;
            // CudaManager& operator=(const CudaManager&) = delete;

            /**
             * @brief Allow shallow copy constructor, but warn
             */
            CudaManager(const CudaManager& other);

            /**
             * @brief Allow shallow assignment operator, but warn
             */
            CudaManager& operator=(const CudaManager& other);


            /**
             * @brief Allow move
             */
            CudaManager(CudaManager&&) = default;
            CudaManager& operator=(CudaManager&&) = default;

            /**
             * @brief Initialize the CUDA environment and select devices.
             * @param device_ids A vector of device IDs to use. If empty, all available devices are used.
             * @throws std::runtime_error if no devices are available or initialization fails.
             */
            void init(const std::vector<int>& device_ids = {});

            /**
             * @brief Get the number of active GPUs.
             * @return The number of active GPUs.
             */
            size_t get_device_count() const;

            /**
             * @brief Get the CUDA stream for a specific device.
             * @param device_id The ID of the device.
             * @return The CUDA stream for the specified device.
             * @throws std::invalid_argument if the device ID is invalid.
             */
            // CUSTREAM get_stream(int device_id) const;

            /**
             * @brief Synchronize all devices.
             */
            void synchronize_all();

            /**
             * @brief Synchronize a specific device.
             * @param device_id The ID of the device to synchronize.
             * @throws std::invalid_argument if the device ID is invalid.
             */
            void synchronize_device(int device_id);

            // /**
            //  * @brief Allocate raw device memory; returns device pointer
            //  */
            // void* allocate(std::size_t size_bytes);

            // /**
            //  * @brief Allocate device memory for custom types
            //  */
            // template <typename T>
            // T* allocate(const T& host_data) {
            //     for (const auto& [device_id, device_manager] : m_device_managers) {
            //         dm->memory().allocate(host_data);
            //     }
            //     std::lock_guard<std::mutex> lock(m_mutex);
            //     T* devptr;
            //     cudaMalloc(&devptr, sizeof(T));
            //     cudaMemcpy(devptr, &host_data, sizeof(T), cudaMemcpyHostToDevice);
            //     m_allocations[devptr] = sizeof(T);
            //     return devptr;
            // };

            /**
             * @brief Create a cuvector for managing device memory on a specific GPU.
             * @tparam T The type of elements in the cuvector.
             * @param device_id The ID of the device.
             * @param size The number of elements to allocate.
             * @return A cuvector of the specified size.
             * @throws std::invalid_argument if the device ID is invalid.
             * @throws std::runtime_error if memory allocation fails.
             */
            // template <template<typename, typename> class VecT, typename T, typename Alloc /* = gpu::CuMAllocator<T> */>
            // VecT<T, Alloc> create_cuvector(int device_id, size_t size);
            // template <typename T>
            // gpu::CUVECTOR_T<T> create_cuvector(int device_id, size_t size);

            /**
             * @brief Copy data from a host vector to a cuvector on a specific GPU.
             * @tparam T The type of elements in the vectors.
             * @param device_id The ID of the device.
             * @param device_vector The cuvector on the device.
             * @param host_vector The host vector containing the data.
             * @throws std::invalid_argument if the device ID is invalid.
             * @throws std::runtime_error if the copy operation fails.
             */
            template <typename T>
            void copy_to_device(int device_id, gpu::CUVECTOR_T<T>& device_vector, const std::vector<T>& host_vector);

            /**
             * @brief Copy data from a cuvector on a specific GPU to a host vector.
             * @tparam T The type of elements in the vectors.
             * @param device_id The ID of the device.
             * @param host_vector The host vector to receive the data.
             * @param device_vector The cuvector on the device.
             * @throws std::invalid_argument if the device ID is invalid.
             * @throws std::runtime_error if the copy operation fails.
             */
            template <typename T>
            void copy_to_host(int device_id, std::vector<T>& host_vector, const gpu::CUVECTOR_T<T>& device_vector);

            /**
             * @brief Get a human-readable description of all active devices.
             * @return A vector of strings describing the active devices.
             */
            std::vector<std::string> get_active_device_descriptions() const;

            /**
             * @brief Automatically select the best CUDA device based on properties.
             * @return The ID of the selected device.
             * @throws std::runtime_error if no suitable device is found.
             */
            int select_best_device() const;

            /**
             * @brief Publish any currently GPU-only ("dirty") fields in
             * `fields` back to the CPU-authoritative
             * configuration::Configuration, BEFORE an algorithm that
             * might read/write them on the CPU side runs. This is the
             * half of the freshness tracker that a plain
             * invalidate-after-the-fact can't cover: if a CPU-side
             * algorithm (e.g. a thermostat) runs while a field is
             * GPU-only-dirty (Leap_Frog_Velocity<gpuBackend>'s velocity
             * write, never round-tripped to CPU), it would otherwise
             * silently read/write the stale pre-write CPU value.
             * Called centrally by Algorithm_Sequence::run() BEFORE
             * every algorithm's apply(), using that algorithm's
             * gpu_mirror_touches() -- default MIRROR_ALL, so every
             * ordinary (CPU-side) algorithm is protected with no
             * changes needed. Declared unconditionally (like
             * invalidate_gpu_mirror() below) since
             * Algorithm_Sequence::run() compiles in CPU-only builds
             * too; the non-CUDA .cc gives it a genuine empty-body
             * no-op (not DISABLED_VOID() -- this is called on *every*
             * algorithm on *every* step regardless of accelerator, by
             * design, so logging a critical message here would fire
             * thousands of times over an ordinary CPU-only run; found
             * the hard way via a real end-to-end run, extended_test/).
             * No-op if `conf` has no mirror yet, or nothing in
             * `fields` is currently dirty.
             */
            void flush_gpu_dirty(configuration::Configuration & conf,
                                  unsigned fields = gpu::MIRROR_ALL);

            /**
             * @brief Publish a CPU-side correction to virial_tensor back
             * into the GPU mirror. Counterpart to flush_gpu_dirty(): that
             * one is mirror-to-CPU (for a CPU-only reader), this one is
             * CPU-to-mirror (for a CPU-only *writer* whose result later
             * GPU-native code still needs to accumulate onto or a later
             * flush must not stomp). Molecular_Virial_Interaction is the
             * only caller today -- it corrects conf.current().virial_
             * tensor in place from atomic to molecular virial; without
             * this, that correction never reaches the mirror, and
             * Pressure_Calculation's own later MIRROR_VIRIAL flush
             * (copy_constraint_data_from_device(), an overwrite not a
             * merge) would silently replace it with the stale raw atomic
             * value again. No-op if `conf` has no mirror yet.
             */
            void publish_cpu_virial(configuration::Configuration & conf);

            /**
             * @brief Clear freshness bits on the Configuration mirror
             * (data-level cache-coherence layer on top of the
             * identity-keyed cache, PLAN.md §3.2 follow-up): declares
             * that `fields` may have been written directly on the CPU
             * side (by the algorithm that just ran) and can no longer
             * be trusted resident on the GPU mirror without a resync.
             * Called centrally by Algorithm_Sequence::run() after
             * every algorithm's apply(), using that algorithm's
             * gpu_mirror_touches() -- default MIRROR_ALL, so every
             * ordinary (CPU-side) algorithm needs no changes to be
             * handled correctly. Only touches gpu_fresh_fields, not
             * gpu_dirty_fields -- flush_gpu_dirty() above already
             * handled publishing anything dirty before this algorithm
             * ran, so by construction nothing it touched should still
             * be dirty by the time this runs. Declared unconditionally
             * (unlike configuration_view()/mark_gpu_dirty() below)
             * because Algorithm_Sequence::run() compiles in CPU-only
             * builds too; the non-CUDA .cc gives it a genuine
             * empty-body no-op (see flush_gpu_dirty()'s comment above
             * for why not DISABLED_VOID()). No-op if `conf` has no
             * mirror yet (nothing to invalidate).
             */
            void invalidate_gpu_mirror(configuration::Configuration & conf,
                                        unsigned fields = gpu::MIRROR_ALL);

            /**
             * @brief Discard `fields`' accumulated producer events
             * (gpu::Configuration::field_producer_events) WITHOUT
             * touching freshness/dirty bits -- unlike invalidate_gpu_
             * mirror(), which also clears gpu_fresh_fields (forcing the
             * next configuration_view() request to think the field is
             * missing and re-upload it from CPU). For a field with no
             * other "starts fresh" point in its lifecycle (e.g.
             * gpu::MIRROR_LATTICE_SHIFT, deliberately excluded from
             * MIRROR_ALL -- see mirror_fields.h -- so nothing generic
             * ever invalidates it), the writer must call this itself
             * once per write, right after its own configuration_view()
             * call has already resynced/waited on whatever was
             * pending, or field_producer_events grows by one entry per
             * step forever -- every future wait_on_fresh_producers()
             * call then has to cudaStreamWaitEvent() the entire ever-
             * growing list, an O(steps^2) cost invisible in short runs
             * (found via a 10000-step benchmark going from 0.8s to
             * 17.5s on this one algorithm's own TIMING line before this
             * method existed). Declared unconditionally like invalidate_
             * gpu_mirror() (empty-body no-op in the non-CUDA .cc); no-op
             * if `conf` has no mirror yet.
             */
            void clear_stale_producer_events(configuration::Configuration & conf,
                                              unsigned fields);

#ifdef USE_CUDA
            /**
             * @brief Identity-keyed GPU mirror cache for topology::Topology
             * (PLAN.md §3.2). Builds the mirror on first call for a given
             * topo.id(); returns the cached one on subsequent calls unless
             * `force_resync` is set -- topology data is static for a
             * normal run, so this is rare (lambda/perturbation-topology
             * updates are the exception). A cache entry that belonged to a
             * different, no-longer-live Topology which happened to be
             * destroyed and have its heap address reused is never a risk
             * here: the cache is keyed on `id()`, a process-wide token
             * that's never reused, not on the object's address -- an
             * unrelated object at the same address has a different id and
             * is correctly treated as a cache miss, not a stale hit.
             */
            gpu::Topology::View topology_view(const topology::Topology & topo,
                                               bool force_resync = false);

            /**
             * @brief Identity-keyed GPU mirror cache for
             * configuration::Configuration (PLAN.md §3.2), now with
             * data-level freshness tracking instead of hand-picked sync
             * booleans: `read_fields` (gpu::MirrorField bits) are the
             * fields the caller needs valid on the GPU side. Any
             * requested field not already marked fresh
             * (gpu::Configuration::gpu_fresh_fields) gets resynced from
             * CPU -- pos+vel via the cheap copy_pos_vel_to_device() path
             * if only POS/VEL were missing, or a full copy_to_device()
             * if FORCE/BOX were requested and missing (no per-field
             * upload routine exists for those, so this falls back to
             * the coarse-grained copy). Freshly synced fields are
             * marked fresh afterwards. Builds the mirror (full sync) on
             * first call for a given conf.id().
             *
             * `stream`, if non-null: fields that are already fresh (no
             * CPU resync needed) but were produced by a GPU kernel on a
             * *different* stream (recorded by mark_gpu_dirty() below)
             * get a cudaStreamWaitEvent(stream, ...) inserted instead of
             * being silently trusted -- this is what makes cross-
             * algorithm ordering safe without a CPU-blocking sync when
             * two GPU-native algorithms run on separate streams (e.g.
             * NonBonded writing MIRROR_FORCE on its own stream, then
             * Leap_Frog_Velocity reading it on its own). Omit (or pass
             * nullptr/0) for legacy/CPU-side callers -- behaves exactly
             * as before, no event bookkeeping.
             */
            gpu::Configuration::View configuration_view(configuration::Configuration & conf,
                                                          unsigned read_fields,
                                                          cudaStream_t stream = 0);

            /**
             * @brief The caller just wrote `fields` into the
             * Configuration mirror via a kernel and is vouching for
             * them being correct/fresh -- no CPU round trip needed.
             * Sets both gpu_fresh_fields (trustworthy, don't
             * re-download) and gpu_dirty_fields (CPU hasn't seen this
             * value yet -- flush_gpu_dirty() must publish it before
             * any CPU-side algorithm touches it). The writer is
             * responsible for also overriding gpu_mirror_touches() so
             * Algorithm_Sequence::run() doesn't immediately undo the
             * fresh bit via invalidate_gpu_mirror() right after it
             * returns.
             *
             * `stream`, if non-null: records one cudaEvent_t (ordering
             * only, cudaEventDisableTiming) marking "everything this
             * call's kernels enqueued on `stream` up to now is done,"
             * and appends it to each written field's producer-event
             * list (gpu::Configuration::field_producer_events) instead
             * of clearing prior entries -- multiple GPU writers can
             * touch the same field across a step (e.g. every bonded
             * term plus NonBonded all write MIRROR_FORCE, each on its
             * own stream) and a later reader must wait on all of them.
             * The list is cleared only at well-defined "this field
             * starts fresh" points: zero_mirror_force(),
             * invalidate_gpu_mirror(), and any CPU->GPU resync inside
             * configuration_view() -- never accumulates across steps.
             */
            void mark_gpu_dirty(configuration::Configuration & conf, unsigned fields,
                                 cudaStream_t stream = 0);

            /**
             * @brief Publish the GPU mirror's current positions/velocities
             * back to the CPU-authoritative configuration::Configuration.
             * For GPU-native integrators (Leap_Frog_*<gpuBackend>) that
             * leave their result resident on the GPU mirror across
             * multiple algorithms and only need one sync-back at the very
             * end, rather than after every kernel. No-op if `conf` has
             * never been mirrored (nothing to publish). Clears
             * gpu_dirty_fields for POS/VEL (CPU has now seen them) but
             * leaves gpu_fresh_fields untouched -- the GPU copy is still
             * trustworthy (it now matches CPU exactly).
             */
            void sync_configuration_from_device(configuration::Configuration & conf);
#endif

            /**
             * @brief Swap the GPU mirror's own current/old halves
             * (gpu::Configuration::exchange_state(), an O(1) struct/
             * pointer swap -- no data movement), at the exact same
             * logical point conf.exchange_state() (CPU) is called: both
             * backends of Leap_Frog_Velocity::apply(), nowhere else.
             * Without this, the mirror's `current`/`old` labels never
             * track which physical state each currently holds, since
             * they're never swapped -- only the CPU-authoritative
             * configuration::Configuration was. Declared unconditionally
             * (like zero_mirror_force()) since Leap_Frog_Velocity<
             * cpuBackend>::apply() compiles in CPU-only builds too;
             * genuine no-op there. No-op if `conf` has no mirror yet.
             */
            void exchange_mirror_state(configuration::Configuration & conf);

            /**
             * @brief Zero the GPU mirror's current().force (and
             * current().virial_tensor), matching Forcefield::
             * calculate_interactions()'s own CPU-side `conf.current().
             * force = 0.0` zero, once per step, before any force-
             * computing Interaction (NonBonded, bonded terms) runs.
             * Every such Interaction then atomicAdd's into the mirror's
             * shared force/virial buffers instead of keeping a private
             * scratch buffer and syncing/reading it back individually --
             * same "zero once centrally, then only ever accumulate"
             * discipline as the constraint-error-flags buffer. Declared
             * unconditionally (like flush_gpu_dirty()) since Forcefield::
             * calculate_interactions() compiles in CPU-only builds too;
             * genuine no-op there. No-op if `conf` has no mirror yet
             * (nothing to zero -- the first real force-computing
             * Interaction's own configuration_view() call builds it,
             * already zeroed via cudaMalloc/resize).
             */
            void zero_mirror_force(configuration::Configuration & conf);

            /**
             * @brief Allocate/resize the GPU mirror's energy_* buffers
             * (gpu::Configuration::resize_energy_groups()) to
             * num_groups doubles each. Called once by each GPU-native
             * bonded/special term's calculate_interactions() (guarded
             * by its own "already resized" flag, since this must run
             * after configuration_view() has guaranteed the mirror
             * exists -- unlike zero_mirror_force()/mark_gpu_dirty(),
             * this can't be a safe no-op on a missing mirror the first
             * time it's needed). Idempotent: every caller passes the
             * same num_energy_groups (fixed for the whole run), so only
             * the first call actually allocates. No-op if `conf` has no
             * mirror yet (should not happen given the ordering above;
             * declared unconditionally like the other mirror methods
             * for CPU-only-build compilation, genuine no-op there).
             */
            void ensure_energy_groups(configuration::Configuration & conf, unsigned num_groups);

            /**
             * @brief Raw device pointers to the mirror's energy_*
             * buffers, for a GPU-native bonded/special term's own
             * launch_accumulate_energy() call. Must be called after
             * ensure_energy_groups() (or another configuration_view()-
             * using call) has already guaranteed the mirror exists;
             * returns all-null EnergyMirrorPtrs otherwise.
             */
            gpu::EnergyMirrorPtrs energy_mirror_ptrs(configuration::Configuration & conf);

            /**
             * @brief Zero the GPU mirror's energy_* buffers, matching
             * zero_mirror_force()'s convention -- called once per step,
             * before any GPU-native bonded/special term's
             * calculate_interactions() runs. No-op if `conf` has no
             * mirror yet or its energy_* buffers aren't allocated yet
             * (nothing to zero).
             */
            void zero_mirror_energy(configuration::Configuration & conf);

            /**
             * @brief Zero the deferred constraint-error-flags buffer
             * (gpu/constraint_error_slots.h). Called once per step, at
             * the top of Algorithm_Sequence::run(), before any
             * constraint algorithm's apply() runs. Declared
             * unconditionally (like flush_gpu_dirty()) since
             * Algorithm_Sequence::run() compiles in CPU-only builds
             * too; genuine no-op there.
             */
            void zero_constraint_error_flags();

            /**
             * @brief The one sync + readback for every constraint
             * algorithm's deferred status, called once at the very end
             * of Algorithm_Sequence::run() instead of each algorithm
             * checking (and syncing on) its own flag immediately.
             * Returns true if any *fatal* slot (gpu::
             * constraint_error_slot_is_fatal()) is nonzero; always
             * populates `out_codes` (size gpu::
             * NUM_CONSTRAINT_ERROR_SLOTS) with every slot's raw value,
             * fatal or not, so the caller can also report non-fatal
             * diagnostics (LINCS's rotation counter). No-op (returns
             * false, `out_codes` left empty) if no constraint algorithm
             * has requested the buffer yet this run (CPU-only builds,
             * or a run using no GPU constraints).
             */
            bool check_constraint_error_flags(std::vector<int> & out_codes);

#ifdef USE_CUDA
            /**
             * @brief Raw device pointer to slot `slot`'s int (gpu::
             * ConstraintErrorSlot) in the shared deferred-error-flags
             * buffer -- pass directly as a kernel's existing
             * `error_flag`/`rotation_count` parameter instead of a
             * private per-algorithm buffer. Lazily allocates the
             * buffer (size gpu::NUM_CONSTRAINT_ERROR_SLOTS) on first
             * call.
             */
            int * constraint_error_flag_slot(unsigned slot);
#endif

        private:
            /**
             * @brief Validate a device ID.
             * @param device_id The ID of the device to validate.
             * @throws std::invalid_argument if the device ID is invalid.
             */
            void validate_device_id(int device_id) const;
#ifdef USE_CUDA
            std::unordered_map<int, std::shared_ptr<CudaDeviceManager> > m_device_managers; ///< Managers for each active device.

            std::unordered_map<std::size_t, std::unique_ptr<gpu::Topology> > m_topologies;
            std::unordered_map<std::size_t, std::unique_ptr<gpu::Configuration> > m_configurations;

            // 1-entry fast path for the overwhelmingly common single-
            // topology/single-configuration case, avoiding a hashmap
            // lookup on every call. 0 is never a real id (util::
            // next_identity_token() starts at 1), so it's a safe "nothing
            // cached yet" sentinel.
            std::size_t m_last_topo_id = 0;
            gpu::Topology * m_last_topo_gpu = nullptr;
            std::size_t m_last_conf_id = 0;
            gpu::Configuration * m_last_conf_gpu = nullptr;

            // Deferred constraint-error-flags buffer (gpu/
            // constraint_error_slots.h) -- one persistent, GPU-resident
            // int array shared by every constraint algorithm, lazily
            // allocated on first constraint_error_flag_slot() call, not
            // per-Configuration (errors aren't conf-scoped).
            gpu::cuvector<int> m_constraint_error_flags;
#endif
    };
}
