/**
 * @file cuda_pairlist_algorithm_impl.cu
 * GPU-native implementation of CUDA_Pairlist_Algorithm_Impl.
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "math/boundary_implementation.h"
#include "math/volume.h"
#include "math/gmath.h"
#include "gpu/cuda/math/periodicity.h"

#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/kernels/periodicity.h"
#include "gpu/cuda/interaction/nonbonded/kernels/lj_crf_tiles.h"
#include "gpu/cuda/interaction/nonbonded/kernels/rf_excluded_kernels.h"
#include "gpu/cuda/interaction/nonbonded/kernels/one_four_kernels.h"
#include "gpu/cuda/interaction/nonbonded/kernels/displacement.h"
#include "gpu/cuda/memory/virial_accumulate_kernels.h"
#include "block_pairlist.h"

#include "cuda_pairlist_algorithm_impl.h"

#include "gpu/cuda/utils.h"


#define NUM_THREADS_PER_BLOCK 256

namespace {
    // Adds the (possibly frozen-from-an-earlier-call) long-range force
    // buffer into the GPU-resident mirror's force array -- one atomicAdd
    // triple per atom, since the mirror may also be receiving concurrent
    // atomicAdd writes from bonded terms running on their own streams.
    __global__ void add_force_into_kernel(FPH3_TYPE * __restrict__ dst,
                                           const FPH3_TYPE * __restrict__ src,
                                           unsigned num_atoms) {
        const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= num_atoms) return;
        atomicAdd(&dst[i].x, src[i].x);
        atomicAdd(&dst[i].y, src[i].y);
        atomicAdd(&dst[i].z, src[i].z);
    }

    void launch_add_force_into(FPH3_TYPE * dst, const FPH3_TYPE * src,
                                unsigned num_atoms, cudaStream_t stream) {
        if (num_atoms == 0) return;
        const unsigned blocks = (num_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
        add_force_into_kernel<<<blocks, NUM_THREADS_PER_BLOCK, 0, stream>>>(dst, src, num_atoms);
    }
}

interaction::CUDA_Pairlist_Algorithm_Impl::CUDA_Pairlist_Algorithm_Impl() {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl constructor");
};

interaction::CUDA_Pairlist_Algorithm_Impl::~CUDA_Pairlist_Algorithm_Impl() {
    if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Pairlist_Algorithm_Impl::init(topology::Topology &topo,
    configuration::Configuration &conf,
    simulation::Simulation &sim,
    std::ostream &os,
    bool quiet) {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl::init");

    // Per-atom iac/charge/energy-group are static for a normal run (same
    // assumption as gpu::Topology's exclusion CSR, built once in its
    // constructor) -- build them here rather than in
    // compute_forces_energies(), which runs every step.
    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    m_iac.resize(num_atoms);
    m_charge.resize(num_atoms);
    m_atom_energy_group.resize(num_atoms);
    for (unsigned i = 0; i < num_atoms; ++i) {
        m_iac[i]              = topo.iac(i);
        m_charge[i]           = static_cast<FPL_TYPE>(topo.charge(i));
        m_atom_energy_group[i] = topo.atom_energy_group(i);
    }
    m_num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
    const unsigned num_buckets = m_num_energy_groups * m_num_energy_groups;

    m_num_atoms = num_atoms;
    if (m_stream == 0) cudaStreamCreate(&m_stream);

    m_e_lj.resize(num_buckets);
    m_e_crf.resize(num_buckets);
    m_virial.resize(9);

    m_longrange_force.resize(num_atoms);
    m_e_lj_long.resize(num_buckets);
    m_e_crf_long.resize(num_buckets);
    m_virial_long.resize(9);
    // Zeroed explicitly (not left as whatever cudaMallocManaged handed
    // back) since these are read every single call, but only written on
    // a recompute_long call -- the very first call, before any
    // classification has run, must see zero long-range contribution.
    cudaMemset(m_longrange_force.data(), 0, sizeof(FPH3_TYPE) * num_atoms);
    cudaMemset(m_e_lj_long.data(),  0, sizeof(double) * num_buckets);
    cudaMemset(m_e_crf_long.data(), 0, sizeof(double) * num_buckets);
    cudaMemset(m_virial_long.data(), 0, sizeof(double) * 9);

    m_candidate_ref_pos.resize(num_atoms);
    m_candidates_built = false;

    // Plain-exclusion CSR for the RF-excluded-pairs kernel
    // (gpu::launch_rf_excluded): built from topo.exclusion(i), NOT
    // topo.all_exclusion(i) (which gpu::Topology's excl_ptr/excl_list
    // uses for tile-classification masking) -- all_exclusion is
    // exclusion() UNION one_four_pair() (topology.cc's
    // update_all_exclusion()), and 1-4 pairs get their own separate
    // lj_exception treatment, never the RF correction. Reusing the
    // tile-classification CSR here silently applied the RF correction to
    // 1-4 pairs too -- caught by rf_excluded_gpu.t.cc's very first run
    // (3x too-negative crf energy, garbled per-atom forces).
    const unsigned num_solute_atoms = topo.num_solute_atoms();
    {
        std::vector<int> h_rf_excl_ptr(num_solute_atoms + 1);
        std::vector<int> h_rf_excl_list;
        h_rf_excl_list.reserve(num_solute_atoms * 3);
        for (unsigned i = 0; i < num_solute_atoms; ++i) {
            h_rf_excl_ptr[i] = static_cast<int>(h_rf_excl_list.size());
            for (topology::excl_cont_t::value_type::const_iterator
                     it = topo.exclusion(i).begin(), to = topo.exclusion(i).end();
                 it != to; ++it) {
                h_rf_excl_list.push_back(static_cast<int>(*it));
            }
        }
        h_rf_excl_ptr[num_solute_atoms] = static_cast<int>(h_rf_excl_list.size());

        m_rf_excl_ptr.resize(num_solute_atoms + 1);
        cudaMemcpy(m_rf_excl_ptr.data(), h_rf_excl_ptr.data(),
                   sizeof(int) * (num_solute_atoms + 1), cudaMemcpyHostToDevice);
        m_rf_excl_list.resize(h_rf_excl_list.size() > 0 ? h_rf_excl_list.size() : 1);
        if (!h_rf_excl_list.empty()) {
            cudaMemcpy(m_rf_excl_list.data(), h_rf_excl_list.data(),
                       sizeof(int) * h_rf_excl_list.size(), cudaMemcpyHostToDevice);
        }
    }

    // 1,4-pair CSR for gpu::launch_one_four -- see cuda_pairlist_algorithm_impl.h's
    // m_one_four_ptr/list doc comment. Same construction as m_rf_excl_ptr/list
    // just above, from topo.one_four_pair(i) instead of topo.exclusion(i).
    {
        std::vector<int> h_one_four_ptr(num_solute_atoms + 1);
        std::vector<int> h_one_four_list;
        h_one_four_list.reserve(num_solute_atoms * 3);
        for (unsigned i = 0; i < num_solute_atoms; ++i) {
            h_one_four_ptr[i] = static_cast<int>(h_one_four_list.size());
            for (topology::excl_cont_t::value_type::const_iterator
                     it = topo.one_four_pair(i).begin(), to = topo.one_four_pair(i).end();
                 it != to; ++it) {
                h_one_four_list.push_back(static_cast<int>(*it));
            }
        }
        h_one_four_ptr[num_solute_atoms] = static_cast<int>(h_one_four_list.size());

        m_one_four_ptr.resize(num_solute_atoms + 1);
        cudaMemcpy(m_one_four_ptr.data(), h_one_four_ptr.data(),
                   sizeof(int) * (num_solute_atoms + 1), cudaMemcpyHostToDevice);
        m_one_four_list.resize(h_one_four_list.size() > 0 ? h_one_four_list.size() : 1);
        if (!h_one_four_list.empty()) {
            cudaMemcpy(m_one_four_list.data(), h_one_four_list.data(),
                       sizeof(int) * h_one_four_list.size(), cudaMemcpyHostToDevice);
        }
    }

    return 0;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::set_cutoff(
                    double const cutoff_short,
                    double const cutoff_long)
{
    m_cutoff_long = cutoff_long;
    m_cutoff_short = cutoff_short;
    m_cutoff_short_2 = cutoff_short * cutoff_short;
    m_cutoff_long_2  = cutoff_long * cutoff_long;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::prepare_cog(
                                configuration::Configuration & conf,
                                topology::Topology & topo,
                                simulation::Simulation & sim) {
    set_cutoff(sim.param().pairlist.cutoff_short,
	     sim.param().pairlist.cutoff_long);

    // Position freshness is needed every step regardless of cutoff mode --
    // reorder()/build_candidates()/classify_tiles() (chargegroup- or
    // atomic-cutoff) all read the GPU position mirror on the next
    // rebuild, and atomic_cutoff mode has no other call site that
    // refreshes it (chargegroup-cutoff mode's box-wrap below needs it too).
    // PLAN.md §3.2: fetched from CudaManager's identity-keyed cache, not
    // conf's own (now-removed) m_gpu member.
    const gpu::Configuration::View conf_view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS);

    if (!sim.param().pairlist.atomic_cutoff){
        const size_t num_cg = topo.num_chargegroups();
        // Sized to *all* chargegroups, not just solute: prepare_cog_kernel
        // (periodicity.cu's prepare_chargegroup) writes cg_cog[cg_i] for
        // every cg_i < num_chargegroups (solute cog is a true centre of
        // geometry, solvent gets its first atom's position -- see that
        // function's doc comment), and classify_tiles() below reads
        // m_cg_cog for solvent candidates too, not just solute ones. A
        // num_solute_chargegroups()-sized buffer here is an out-of-bounds
        // write on every call whenever any solvent chargegroups exist --
        // found via compute-sanitizer while investigating KNOWN_ISSUES.md's
        // "CUDA context corruption after runtime atomic_cutoff toggle"
        // entry: this OOB write happens on every prepare_cog() call
        // regardless of atomic_cutoff, but only sometimes corrupts memory
        // that mattered, which is why it looked toggle-specific.
        m_cg_cog.resize(num_cg);
        m_cg_cells.resize(num_cg);
        m_cg_sort_key.resize(num_cg);
        const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);
        SPLIT_BOUNDARY(_prepare_cog, conf, topo, topo_view, conf_view);
    }
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_prepare_cog<B>(
                                    configuration::Configuration & conf,
                                    topology::Topology & topo,
                                    gpu::Topology::View topo_view,
                                    gpu::Configuration::View conf_view) {
    DEBUG(10, "putting chargegroups into box");

    const unsigned num_cg = static_cast<unsigned>(topo.num_chargegroups());
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((num_cg + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);

    // conf_view/topo_view were fetched (with the one real per-step
    // resync) by prepare_cog(), the only caller.
    gpu::Periodicity<B> periodicity(conf.current().box);
    periodicity.set_cell_size(m_cutoff_long);
    gpu::prepare_cog_kernel<<<dimGrid, dimBlock>>>(topo_view,
                                                    conf_view,
                                                    periodicity,
                                                    m_cg_cog.view(),
                                                    m_cg_cells.view(),
                                                    m_cg_sort_key.view());
};

void interaction::CUDA_Pairlist_Algorithm_Impl::reorder(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    const unsigned num_solute_atoms  = static_cast<unsigned>(topo.num_solute_atoms());
    const unsigned num_atoms         = static_cast<unsigned>(topo.num_atoms());
    const unsigned num_solvent_atoms = num_atoms - num_solute_atoms;
    const unsigned num_cg            = static_cast<unsigned>(topo.num_chargegroups());

    const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

    // MIRROR_POS: prepare_cog() already requested (and, in chargegroup-
    // cutoff mode, box-wrapped) it this step, marking it fresh on the
    // GPU mirror -- this resolves to a cache hit, no re-download, so
    // the box-wrapped GPU positions aren't overwritten by the
    // (unwrapped) CPU copy before reorder() ever reads them.
    // Not const: ConfigurationView::current()/old() aren't const-qualified.
    gpu::Configuration::View conf_view =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS);

    m_atom_sort_key.resize(num_atoms);
    if (!sim.param().pairlist.atomic_cutoff) {
        // Every atom inherits its owning chargegroup's Morton cell key
        // (TILE_PAIRLIST_DESIGN.md §3 step 1) -- keeps chargegroup members
        // spatially adjacent after sorting, without a separate per-atom
        // Morton computation.
        dim3 dimBlock(NUM_THREADS_PER_BLOCK);
        dim3 dimGrid((num_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
        gpu::atom_sort_key_kernel<<<dimGrid, dimBlock>>>(
            topo_view.chargegroup, num_cg, num_atoms,
            m_cg_sort_key.data(), m_atom_sort_key.data());
    } else {
        // atomic_cutoff mode has no chargegroup-level cog/cell build to
        // inherit from (prepare_cog() skips it entirely) -- compute each
        // atom's own Morton key directly (TILE_PAIRLIST_DESIGN.md §4.2/
        // step 7).
        SPLIT_BOUNDARY(_atom_sort_key_atomic, conf, num_atoms, conf_view);
    }

    m_solute_atom_order.resize(num_solute_atoms);
    m_solvent_atom_order.resize(num_solvent_atoms);

    // Identity permutation, in global atom numbering (solute first, then
    // solvent -- matching m_atom_sort_key and conf's atom position array).
    // Hand-written kernel, not thrust::sequence -- see sequence_kernel's
    // doc comment (block_pairlist.h) for why.
    if (num_solute_atoms > 0) {
        const unsigned blocks = (num_solute_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
        gpu::sequence_kernel<<<blocks, NUM_THREADS_PER_BLOCK>>>(m_solute_atom_order.data(), num_solute_atoms, 0u);
    }
    if (num_solvent_atoms > 0) {
        const unsigned blocks = (num_solvent_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;
        gpu::sequence_kernel<<<blocks, NUM_THREADS_PER_BLOCK>>>(m_solvent_atom_order.data(), num_solvent_atoms, num_solute_atoms);
    }

    // Sort each group's identity permutation by its (inherited) Morton
    // cell key. m_atom_sort_key is mutated in place -- safe, it's
    // recomputed fresh above every call. Hand-written bitonic sort, not
    // thrust::sort_by_key -- see sequence_kernel's doc comment
    // (block_pairlist.h) for why.
    if (num_solute_atoms > 0) {
        gpu::bitonic_sort_by_key(m_atom_sort_key.data(), m_solute_atom_order.data(),
                                  num_solute_atoms, m_sort_scratch_keys, m_sort_scratch_values);
    }
    if (num_solvent_atoms > 0) {
        gpu::bitonic_sort_by_key(m_atom_sort_key.data() + num_solute_atoms, m_solvent_atom_order.data(),
                                  num_solvent_atoms, m_sort_scratch_keys, m_sort_scratch_values);
    }

    const unsigned num_solute_blocks  = (num_solute_atoms  + gpu::BLOCK_SIZE - 1) / gpu::BLOCK_SIZE;
    const unsigned num_solvent_blocks = (num_solvent_atoms + gpu::BLOCK_SIZE - 1) / gpu::BLOCK_SIZE;

    m_solute_block_center.resize(num_solute_blocks);
    m_solute_block_radius.resize(num_solute_blocks);
    m_solvent_block_center.resize(num_solvent_blocks);
    m_solvent_block_radius.resize(num_solvent_blocks);

    // Bounding spheres come directly from each atom's own position now
    // (box-wrapped in chargegroup-cutoff mode, raw/unwrapped in
    // atomic-cutoff mode -- see _atom_sort_key_atomic's doc comment), not
    // a chargegroup cog -- m_cg_cog is reserved for classify_tiles()'s
    // chargegroup-cutoff distance test.
    const math::CuVArray::View pos = conf_view.current().pos;

    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    if (num_solute_blocks > 0) {
        dim3 dimGrid((num_solute_blocks + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
        gpu::compute_block_bounds_kernel<<<dimGrid, dimBlock>>>(
            m_solute_atom_order.data(), num_solute_atoms, num_solute_blocks,
            pos, m_solute_block_center.data(), m_solute_block_radius.data());
    }
    if (num_solvent_blocks > 0) {
        dim3 dimGrid((num_solvent_blocks + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
        gpu::compute_block_bounds_kernel<<<dimGrid, dimBlock>>>(
            m_solvent_atom_order.data(), num_solvent_atoms, num_solvent_blocks,
            pos, m_solvent_block_center.data(), m_solvent_block_radius.data());
    }
};

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_atom_sort_key_atomic(
                configuration::Configuration & conf,
                unsigned num_atoms,
                gpu::Configuration::View conf_view)
{
    gpu::Periodicity<B> periodicity(conf.current().box);
    // get_cell() requires set_cell_size() first (needs a real cutoff-sized
    // cell, not the default-constructed zero size).
    periodicity.set_cell_size(m_cutoff_long);

    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((num_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
    gpu::atom_cell_kernel<B><<<dimGrid, dimBlock>>>(
        conf_view, num_atoms, periodicity, m_atom_sort_key.view());
};

void interaction::CUDA_Pairlist_Algorithm_Impl::build_candidates(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    // v1 scope (TILE_PAIRLIST_DESIGN.md §4.1): vacuum + rectangular only.
    // CUDA_Pairlist_Algorithm::init() already hard-errors for other
    // boundary types before this is ever reached; SPLIT_BOUNDARY still
    // needs a concrete case for every math::boundary_enum value, so
    // _build_candidates<triclinic/truncoct> exist but are unreachable.
    SPLIT_BOUNDARY(_build_candidates, conf, topo, sim);
}

namespace {
  /**
   * Rough capacity estimate for a candidate TileVecT, mirroring
   * Nonbonded_Set::init's per-atom pairlist.reserve() density estimate
   * (nonbonded_set.cc) but in terms of blocks. This is a sizing heuristic
   * only, not a correctness mechanism -- the real safeguard is the
   * was_overflowed() check after the kernel runs (see below). Falls back
   * to the exact worst case (every block pair) for vacuum / degenerate
   * volumes, where no density estimate is meaningful, and never exceeds
   * that worst case otherwise either.
   */
  unsigned estimate_candidate_capacity(unsigned num_blocks_a, unsigned num_blocks_b,
                                        double cutoff, double vol, bool self_pairs) {
    const unsigned worst_case = self_pairs
        ? (num_blocks_a * (num_blocks_a + 1)) / 2
        : num_blocks_a * num_blocks_b;
    if (vol <= 0.0 || worst_case == 0) return worst_case;

    const double block_density = num_blocks_b / vol;
    const double sphere_vol = 4.0 / 3.0 * math::Pi * cutoff * cutoff * cutoff;
    // 1.3x matches the CPU estimate's safety margin; the extra 2x accounts
    // for block bounding spheres being a conservative (over-inclusive)
    // test compared to the CPU's exact per-atom distance check.
    double per_block = 1.3 * 2.0 * block_density * sphere_vol;
    if (self_pairs) per_block *= 0.5; // only bi <= bj pairs get pushed
    // ceil(), not a truncating cast: a truncating cast silently rounds
    // e.g. 0.93 down to 0, systematically UNDERestimating capacity by up
    // to 1 per block -- found via the pairlist-equivalence test
    // (TILE_PAIRLIST_DESIGN.md §5/§6): a real, reproducible overflow (not
    // just an off-by-one on paper) that TileVecT correctly detected and
    // reported via was_overflowed(), but the report went to io::messages
    // and nothing had displayed it yet, so it looked like silent data
    // corruption instead of the flagged error it actually was.
    const unsigned capacity = static_cast<unsigned>(std::ceil(per_block * num_blocks_a)) + num_blocks_a;
    return std::min(capacity, worst_case);
  }

  bool candidate_overflowed(gpu::TileVecT<gpu::Interaction_Tile> const & tiles) {
    cudaDeviceSynchronize();
    return tiles.was_overflowed();
  }
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_build_candidates(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    const unsigned num_solute_blocks  = static_cast<unsigned>(m_solute_block_center.size());
    const unsigned num_solvent_blocks = static_cast<unsigned>(m_solvent_block_center.size());

    m_tiles.clear();

    gpu::Periodicity<B> periodicity(conf.current().box);
    const double cutoff_d = m_cutoff_long + sim.param().pairlist.skin;
    const FPL_TYPE cutoff = static_cast<FPL_TYPE>(cutoff_d);
    const double vol = math::volume(conf.current().box, conf.boundary_type);

    dim3 dimBlock2D(16, 16);

    // solute_candidates is shared by two passes below (solute-solute and
    // solute-solvent), both writing through the same atomically-advanced
    // m_size -- reserve their combined capacity once, up front. TileVecT::
    // reserve() deallocates-and-reallocates when growing, which would
    // silently discard the first pass's results if called again in
    // between; call it exactly once here, before either pass runs.
    auto run_solute_candidates = [&]() {
        if (num_solute_blocks > 0) {
            // solute - solute (self-pairs, bi <= bj only)
            dim3 dimGrid((num_solute_blocks + 15) / 16, (num_solute_blocks + 15) / 16);
            gpu::find_block_candidates_kernel<B><<<dimGrid, dimBlock2D>>>(
                num_solute_blocks, num_solute_blocks,
                m_solute_block_center.data(), m_solute_block_radius.data(),
                m_solute_block_center.data(), m_solute_block_radius.data(),
                true, periodicity, cutoff, m_tiles.solute_candidates.view());
        }
        if (num_solute_blocks > 0 && num_solvent_blocks > 0) {
            // solute - solvent (every pair; column side is solvent-order)
            dim3 dimGrid((num_solute_blocks + 15) / 16, (num_solvent_blocks + 15) / 16);
            gpu::find_block_candidates_kernel<B><<<dimGrid, dimBlock2D>>>(
                num_solute_blocks, num_solvent_blocks,
                m_solute_block_center.data(), m_solute_block_radius.data(),
                m_solvent_block_center.data(), m_solvent_block_radius.data(),
                false, periodicity, cutoff, m_tiles.solute_candidates.view());
        }
    };

    if (num_solute_blocks > 0) {
        const unsigned cap_ss = estimate_candidate_capacity(num_solute_blocks, num_solute_blocks, cutoff_d, vol, true);
        const unsigned cap_sv = num_solvent_blocks > 0
            ? estimate_candidate_capacity(num_solute_blocks, num_solvent_blocks, cutoff_d, vol, false)
            : 0;
        m_tiles.solute_candidates.reserve(cap_ss + cap_sv);
        run_solute_candidates();

        // The density estimate above is just a sizing heuristic (see its
        // doc comment) -- a real system (extended_test/ubiquitin) can and
        // does exceed it, since block bounding spheres don't have a
        // uniform radius the way the estimate assumes. Retry once at the
        // true worst-case capacity (every possible block pair -- always
        // mathematically sufficient, no further retry can ever be
        // needed) rather than silently dropping candidates and running
        // with a wrong pairlist. Still cheap in the common (no overflow)
        // case: this branch only runs at all when was_overflowed() is
        // actually set.
        if (candidate_overflowed(m_tiles.solute_candidates)) {
            const unsigned worst_case =
                (num_solute_blocks * (num_solute_blocks + 1)) / 2 +
                num_solute_blocks * num_solvent_blocks;
            m_tiles.solute_candidates.reserve(worst_case);
            run_solute_candidates();
            if (candidate_overflowed(m_tiles.solute_candidates)) {
                // Unreachable in practice -- worst_case is a hard upper
                // bound on how many block pairs can possibly exist.
                // Reaching this means TileVecT itself is broken, not the
                // sizing heuristic.
                io::messages.add(
                    "CUDA pairlist candidate build overflowed solute_candidates "
                    "even at worst-case capacity -- this indicates a bug in "
                    "TileVecT, not the sizing heuristic.",
                    "CUDA_Pairlist_Algorithm", io::message::error);
            }
        }
    }

    if (num_solvent_blocks > 0) {
        auto run_solvent_candidates = [&]() {
            // solvent - solvent (self-pairs, bi <= bj only)
            dim3 dimGrid((num_solvent_blocks + 15) / 16, (num_solvent_blocks + 15) / 16);
            gpu::find_block_candidates_kernel<B><<<dimGrid, dimBlock2D>>>(
                num_solvent_blocks, num_solvent_blocks,
                m_solvent_block_center.data(), m_solvent_block_radius.data(),
                m_solvent_block_center.data(), m_solvent_block_radius.data(),
                true, periodicity, cutoff, m_tiles.solvent_candidates.view());
        };

        m_tiles.solvent_candidates.reserve(
            estimate_candidate_capacity(num_solvent_blocks, num_solvent_blocks, cutoff_d, vol, true));
        run_solvent_candidates();

        if (candidate_overflowed(m_tiles.solvent_candidates)) {
            const unsigned worst_case = (num_solvent_blocks * (num_solvent_blocks + 1)) / 2;
            m_tiles.solvent_candidates.reserve(worst_case);
            run_solvent_candidates();
            if (candidate_overflowed(m_tiles.solvent_candidates)) {
                io::messages.add(
                    "CUDA pairlist candidate build overflowed solvent_candidates "
                    "even at worst-case capacity -- this indicates a bug in "
                    "TileVecT, not the sizing heuristic.",
                    "CUDA_Pairlist_Algorithm", io::message::error);
            }
        }
    }
};

bool interaction::CUDA_Pairlist_Algorithm_Impl::needs_candidate_rebuild(
                configuration::Configuration & conf,
                simulation::Simulation & sim)
{
    if (!m_candidates_built) return true;

    // Explicit early-return, not left as an emergent property of the
    // comparison below -- see this method's doc comment (header).
    const double skin = sim.param().pairlist.skin;
    if (skin == 0.0) return true;

    const unsigned num_atoms = static_cast<unsigned>(m_candidate_ref_pos.size());
    // sync_pos_vel = false: prepare() already did this step's one real
    // resync; called from update(), always after prepare().
    const math::CuVArray::View current_pos =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;
    const FPL_TYPE max_disp = gpu::launch_max_displacement(
        current_pos, m_candidate_ref_pos.view(), num_atoms,
        conf.boundary_type, conf.current().box, m_displacement_partial);

    return static_cast<double>(2.0 * max_disp) >= skin;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::rebuild_candidates(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    reorder(conf, topo, sim);
    build_candidates(conf, topo, sim);

    // Snapshot positions *after* reorder()/build_candidates() ran (they
    // don't mutate positions themselves, but this keeps the "as of last
    // rebuild" semantics unambiguous regardless).
    const unsigned num_atoms = static_cast<unsigned>(m_candidate_ref_pos.size());
    const math::CuVArray::View current_pos =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;
    for (unsigned i = 0; i < num_atoms; ++i) {
        m_candidate_ref_pos[i] = current_pos(i);
    }
    m_candidates_built = true;
    ++m_candidate_rebuild_count;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::classify_tiles(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    // Same v1 boundary-scope note as _build_candidates: init() already
    // hard-errors for triclinic/truncoct, but SPLIT_BOUNDARY still needs
    // a concrete case for every math::boundary_enum value.
    SPLIT_BOUNDARY(_classify_tiles, conf, topo, sim);
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_classify_tiles(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    m_tiles.solute_short.clear();
    m_tiles.solute_long.clear();
    m_tiles.solvent_short.clear();
    m_tiles.solvent_long.clear();

    // At most one output tile per input candidate tile per bucket -- an
    // exact bound, unlike the candidate search's own density estimate
    // (which can't know the true block-pair count ahead of running it).
    const unsigned num_solute_candidates  = static_cast<unsigned>(m_tiles.solute_candidates.size());
    const unsigned num_solvent_candidates = static_cast<unsigned>(m_tiles.solvent_candidates.size());
    if (num_solute_candidates > 0) {
        m_tiles.solute_short.reserve(num_solute_candidates);
        m_tiles.solute_long.reserve(num_solute_candidates);
    }
    if (num_solvent_candidates > 0) {
        m_tiles.solvent_short.reserve(num_solvent_candidates);
        m_tiles.solvent_long.reserve(num_solvent_candidates);
    }

    gpu::Periodicity<B> periodicity(conf.current().box);
    const FPL_TYPE cutoff_short2 = static_cast<FPL_TYPE>(m_cutoff_short_2);
    const FPL_TYPE cutoff_long2  = static_cast<FPL_TYPE>(m_cutoff_long_2);
    const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);
    // sync_pos_vel = false: prepare() already did this step's one real
    // resync; called from update(), always after prepare().
    const math::CuVArray::View pos =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;
    const bool atomic_cutoff = sim.param().pairlist.atomic_cutoff;

    dim3 dimBlock2D(gpu::BLOCK_SIZE, gpu::BLOCK_SIZE);

    if (num_solute_candidates > 0) {
        if (atomic_cutoff) {
            gpu::classify_tiles_kernel<true, B><<<num_solute_candidates, dimBlock2D>>>(
                m_tiles.solute_candidates.view(),
                m_solute_atom_order.data(), static_cast<unsigned>(m_solute_atom_order.size()),
                m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
                m_cg_cog.view(), pos, topo_view, periodicity,
                cutoff_short2, cutoff_long2,
                m_tiles.solute_short.view(), m_tiles.solute_long.view());
        } else {
            gpu::classify_tiles_kernel<false, B><<<num_solute_candidates, dimBlock2D>>>(
                m_tiles.solute_candidates.view(),
                m_solute_atom_order.data(), static_cast<unsigned>(m_solute_atom_order.size()),
                m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
                m_cg_cog.view(), pos, topo_view, periodicity,
                cutoff_short2, cutoff_long2,
                m_tiles.solute_short.view(), m_tiles.solute_long.view());
        }
    }

    if (num_solvent_candidates > 0) {
        if (atomic_cutoff) {
            gpu::classify_tiles_kernel<true, B><<<num_solvent_candidates, dimBlock2D>>>(
                m_tiles.solvent_candidates.view(),
                m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
                nullptr, 0u,
                m_cg_cog.view(), pos, topo_view, periodicity,
                cutoff_short2, cutoff_long2,
                m_tiles.solvent_short.view(), m_tiles.solvent_long.view());
        } else {
            gpu::classify_tiles_kernel<false, B><<<num_solvent_candidates, dimBlock2D>>>(
                m_tiles.solvent_candidates.view(),
                m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
                nullptr, 0u,
                m_cg_cog.view(), pos, topo_view, periodicity,
                cutoff_short2, cutoff_long2,
                m_tiles.solvent_short.view(), m_tiles.solvent_long.view());
        }
    }

    cudaDeviceSynchronize();
    if (m_tiles.solute_short.was_overflowed()  || m_tiles.solute_long.was_overflowed() ||
        m_tiles.solvent_short.was_overflowed() || m_tiles.solvent_long.was_overflowed()) {
        io::messages.add(
          "CUDA pairlist classification overflowed a short/long tile "
          "capacity -- this should be impossible (reserved capacity was "
          "the exact candidate count); a real bug, not a sizing-heuristic "
          "issue.",
          "CUDA_Pairlist_Algorithm", io::message::error);
    }
};

namespace {
  /**
   * Unpack every tile in `tiles` into (atom_i, atom_j) pairs (i < j) and
   * append j to out[i], using the same row/col atom-order convention as
   * classify_tiles_kernel: col_other_order is only consulted when a
   * tile's col_from_b flag is set (the solute-solvent case); pass
   * nullptr for it when unpacking solvent_short/solvent_long, where that
   * flag is never set.
   */
  void unpack_tiles_into(const gpu::TileVecT<gpu::Interaction_Tile> & tiles,
                          const gpu::cuvector<unsigned> & row_order,
                          const gpu::cuvector<unsigned> * col_other_order,
                          interaction::Pairlist & out) {
    const unsigned n = static_cast<unsigned>(tiles.size());
    for (unsigned t = 0; t < n; ++t) {
        const gpu::Interaction_Tile & tile = tiles[t];
        unsigned row_block, col_block;
        bool col_from_b;
        gpu::unpack_block_index(tile.index, row_block, col_block, col_from_b);
        const gpu::cuvector<unsigned> & col_order =
            (col_from_b && col_other_order) ? *col_other_order : row_order;

        for (unsigned r = 0; r < gpu::BLOCK_SIZE; ++r) {
            const unsigned bits = tile.mask[r];
            if (!bits) continue;
            const unsigned row_idx = row_block * gpu::BLOCK_SIZE + r;
            if (row_idx >= row_order.size()) continue;
            const unsigned a1 = row_order[row_idx];
            for (unsigned c = 0; c < gpu::BLOCK_SIZE; ++c) {
                if (!(bits & (1u << c))) continue;
                const unsigned col_idx = col_block * gpu::BLOCK_SIZE + c;
                if (col_idx >= col_order.size()) continue;
                const unsigned a2 = col_order[col_idx];
                const unsigned i = a1 < a2 ? a1 : a2;
                const unsigned j = a1 < a2 ? a2 : a1;
                out[i].push_back(j);
            }
        }
    }
  }
}

interaction::PairlistContainer interaction::CUDA_Pairlist_Algorithm_Impl::to_pairlist_container(
                topology::Topology & topo) const
{
    interaction::PairlistContainer result;
    result.resize(static_cast<unsigned>(topo.num_atoms()));

    // Tiles/order arrays live in unified memory, written by kernels that
    // classify_tiles() already synchronized after -- safe to read here.
    unpack_tiles_into(m_tiles.solute_short,  m_solute_atom_order, &m_solvent_atom_order, result.solute_short);
    unpack_tiles_into(m_tiles.solute_long,   m_solute_atom_order, &m_solvent_atom_order, result.solute_long);
    unpack_tiles_into(m_tiles.solvent_short, m_solvent_atom_order, nullptr, result.solvent_short);
    unpack_tiles_into(m_tiles.solvent_long,  m_solvent_atom_order, nullptr, result.solvent_long);

    return result;
}

void interaction::CUDA_Pairlist_Algorithm_Impl::compute_forces_energies(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim,
                gpu::LJParamView lj,
                gpu::NbSimParams nb,
                bool recompute_long,
                bool rf_excluded)
{
    const unsigned num_atoms   = m_num_atoms;
    const unsigned num_buckets = m_num_energy_groups * m_num_energy_groups;

    // IEEE-754 zero is the all-zero bit pattern -- cudaMemsetAsync is
    // safe and matches the existing convention (TileVecT::clear()
    // zeroes its tile array the same way). Short-range force is no
    // longer a private buffer -- it's written directly into the
    // GPU-resident mirror (zeroed once per step by Forcefield::
    // calculate_interactions()'s sim.cuda().zero_mirror_force(), not
    // here).
    cudaMemsetAsync(m_e_lj.data(),  0, sizeof(double) * num_buckets, m_stream);
    cudaMemsetAsync(m_e_crf.data(), 0, sizeof(double) * num_buckets, m_stream);
    cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

    // sync_pos_vel = false: prepare() already did this step's one real
    // resync; called from calculate_interactions(), always after
    // prepare(). Only MIRROR_POS is requested as a read field (never
    // MIRROR_FORCE) -- that's what keeps this from triggering a coarse
    // resync that would clobber force bonded terms may have already
    // accumulated into the mirror this step.
    gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, m_stream);
    const math::CuVArray::View pos = view.current().pos;
    FPH3_TYPE * const mirror_force = view.current().force.data();
    const math::boundary_enum boundary = conf.boundary_type;
    const math::Box box = conf.current().box;

    const unsigned num_solute_atoms  = static_cast<unsigned>(m_solute_atom_order.size());
    const unsigned num_solvent_atoms = static_cast<unsigned>(m_solvent_atom_order.size());

    // Short-range: real GROMOS twin-range recomputes this every step,
    // regardless of recompute_long -- same row_order/col_other_order
    // convention as classify_tiles_kernel (block_pairlist.h):
    // solute_short/solute_long tiles may reference either solute or
    // solvent atoms on their column side (col_from_b), solvent_short/
    // solvent_long tiles never do.
    gpu::launch_lj_crf_tiles(
        m_tiles.solute_short.view(), m_solute_atom_order.data(), num_solute_atoms,
        m_solvent_atom_order.data(), num_solvent_atoms,
        pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
        mirror_force, m_e_lj.data(), m_e_crf.data(), m_virial.data(), m_stream);

    gpu::launch_lj_crf_tiles(
        m_tiles.solvent_short.view(), m_solvent_atom_order.data(), num_solvent_atoms,
        nullptr, 0u,
        pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
        mirror_force, m_e_lj.data(), m_e_crf.data(), m_virial.data(), m_stream);

    // RF for excluded pairs (param.nonbonded.rf_excluded): independent of
    // the pairlist/twin-range cadence entirely -- walks the exclusion
    // list and solvent chargegroups directly, same as
    // nonbonded_set.cc's unconditional RF_excluded_outerloop call.
    // Accumulates into the same m_force/m_e_crf/m_virial buffers as the
    // short-range tile kernel above, so it's covered by the same
    // conf.current() add-in below.
    if (rf_excluded) {
        const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);
        gpu::launch_rf_excluded(
            m_rf_excl_ptr.data(), m_rf_excl_list.data(), num_solute_atoms,
            topo_view.chargegroup, topo_view.num_solute_chargegroups, topo_view.num_chargegroups,
            pos, m_charge.data(), m_atom_energy_group.data(), nb, boundary, box,
            mirror_force, m_e_crf.data(), m_virial.data(), m_stream);
    }

    // 1,4-pair ("LJ exception") interactions: unconditional, matching
    // nonbonded_set.cc's unconditional one_four_outerloop call -- these
    // pairs are excluded from the tile kernels above (all_exclusion =
    // exclusion UNION one_four_pair) and only ever computed here, with
    // their own scaled cs6/cs12 LJ parameters and Coulomb-scaled CRF.
    // Accumulates into the same m_force/m_e_lj/m_e_crf/m_virial buffers.
    gpu::launch_one_four(
        m_one_four_ptr.data(), m_one_four_list.data(), num_solute_atoms,
        pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb,
        nb.coulomb_scaling, boundary, box,
        mirror_force, m_e_lj.data(), m_e_crf.data(), m_virial.data(), m_stream);

    // Long-range: only recomputed on a classification/rebuild step
    // (recompute_long == true) -- otherwise m_longrange_force/
    // m_e_lj_long/m_e_crf_long are left untouched, holding whatever was
    // last computed here, matching nonbonded_set.cc's frozen
    // m_longrange_storage exactly. The zero must stay inside this
    // branch: zeroing unconditionally would wipe the frozen values a
    // non-rebuild step is supposed to reuse.
    if (recompute_long) {
        cudaMemsetAsync(m_longrange_force.data(), 0, sizeof(FPH3_TYPE) * num_atoms, m_stream);
        cudaMemsetAsync(m_e_lj_long.data(),  0, sizeof(double) * num_buckets, m_stream);
        cudaMemsetAsync(m_e_crf_long.data(), 0, sizeof(double) * num_buckets, m_stream);
        cudaMemsetAsync(m_virial_long.data(), 0, sizeof(double) * 9, m_stream);

        gpu::launch_lj_crf_tiles(
            m_tiles.solute_long.view(), m_solute_atom_order.data(), num_solute_atoms,
            m_solvent_atom_order.data(), num_solvent_atoms,
            pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
            m_longrange_force.data(), m_e_lj_long.data(), m_e_crf_long.data(), m_virial_long.data(), m_stream);

        gpu::launch_lj_crf_tiles(
            m_tiles.solvent_long.view(), m_solvent_atom_order.data(), num_solvent_atoms,
            nullptr, 0u,
            pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
            m_longrange_force.data(), m_e_lj_long.data(), m_e_crf_long.data(), m_virial_long.data(), m_stream);
    }

    // Long-range force is a private, persistent buffer (must freeze
    // across non-recompute_long steps, unlike short-range which is
    // recomputed every call) -- add its current value into the mirror
    // every single step, regardless of recompute_long, since the mirror
    // itself is zeroed fresh every step by Forcefield and would
    // otherwise lose the frozen contribution entirely on a non-rebuild
    // step. Same accumulation nonbonded_set.cc's m_storage.force +=
    // m_longrange_storage.force does, just GPU-resident.
    launch_add_force_into(mirror_force, m_longrange_force.data(), num_atoms, m_stream);

    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE, m_stream);

    // Energy/virial are small, private, double-precision buffers (not
    // part of the mirror) so still need a sync to read back, but only on
    // this algorithm's own stream, not a device-wide barrier.
    cudaStreamSynchronize(m_stream);

    // Direct per-[gi][gj] accumulation, same style as
    // nonbonded_innerloop.cc's CPU inner loop -- no scalar out-params.
    for (unsigned gi = 0; gi < m_num_energy_groups; ++gi) {
        for (unsigned gj = 0; gj < m_num_energy_groups; ++gj) {
            const unsigned k = gi * m_num_energy_groups + gj;
            conf.current().energies.lj_energy[gi][gj]  += m_e_lj[k]  + m_e_lj_long[k];
            conf.current().energies.crf_energy[gi][gj] += m_e_crf[k] + m_e_crf_long[k];
        }
    }

    // Atomic virial, exact CPU formula (nonbonded_innerloop.cc):
    // virial_tensor(b, a) += r(b) * force(a). Not energy-group-bucketed.
    // Whether this is *used* downstream (plain atomic virial, or further
    // corrected to molecular virial by Molecular_Virial_Interaction,
    // generic across accelerators) is entirely Forcefield's decision.
    // Published into the shared mirror's virial_tensor via atomicAdd
    // (GPU-resident, no CPU round trip) -- every bonded term/NonBonded/
    // active constraint algorithm contributes to the same global
    // accumulator each step, zeroed once per step by CudaManager::
    // zero_mirror_force(). Two calls (short-range + frozen long-range),
    // not a pre-summed one: atomicAdd is commutative/associative, and
    // this avoids an extra elementwise host-or-device combine step.
    gpu::launch_accumulate_virial9(
        reinterpret_cast<FPH_TYPE*>(view.current().virial_tensor), m_virial.data(), m_stream);
    gpu::launch_accumulate_virial9(
        reinterpret_cast<FPH_TYPE*>(view.current().virial_tensor), m_virial_long.data(), m_stream);
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VIRIAL, m_stream);
}
