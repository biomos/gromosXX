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
#include "gpu/cuda/interaction/nonbonded/kernels/displacement.h"
#include "block_pairlist.h"

#include "cuda_pairlist_algorithm_impl.h"

#include "gpu/cuda/utils.h"

#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#define NUM_THREADS_PER_BLOCK 256

interaction::CUDA_Pairlist_Algorithm_Impl::CUDA_Pairlist_Algorithm_Impl() {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl constructor");
};

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

    m_force.resize(num_atoms);
    m_e_lj.resize(num_buckets);
    m_e_crf.resize(num_buckets);

    m_longrange_force.resize(num_atoms);
    m_e_lj_long.resize(num_buckets);
    m_e_crf_long.resize(num_buckets);
    // Zeroed explicitly (not left as whatever cudaMallocManaged handed
    // back) since these are read every single call, but only written on
    // a recompute_long call -- the very first call, before any
    // classification has run, must see zero long-range contribution.
    cudaMemset(m_longrange_force.data(), 0, sizeof(FPL3_TYPE) * num_atoms);
    cudaMemset(m_e_lj_long.data(),  0, sizeof(double) * num_buckets);
    cudaMemset(m_e_crf_long.data(), 0, sizeof(double) * num_buckets);

    m_candidate_ref_pos.resize(num_atoms);
    m_candidates_built = false;

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
        const size_t num_solute_cg = topo.num_solute_chargegroups();
        const size_t num_cg = topo.num_chargegroups();
        m_cg_cog.resize(num_solute_cg);
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
    if (num_solute_atoms > 0) {
        thrust::sequence(thrust::device, m_solute_atom_order.data(), m_solute_atom_order.data() + num_solute_atoms, 0u);
    }
    if (num_solvent_atoms > 0) {
        thrust::sequence(thrust::device, m_solvent_atom_order.data(), m_solvent_atom_order.data() + num_solvent_atoms, num_solute_atoms);
    }

    // Sort each group's identity permutation by its (inherited) Morton
    // cell key. m_atom_sort_key is mutated in place -- safe, it's
    // recomputed fresh above every call.
    if (num_solute_atoms > 0) {
        thrust::sort_by_key(thrust::device,
            m_atom_sort_key.data(), m_atom_sort_key.data() + num_solute_atoms,
            m_solute_atom_order.data());
    }
    if (num_solvent_atoms > 0) {
        thrust::sort_by_key(thrust::device,
            m_atom_sort_key.data() + num_solute_atoms, m_atom_sort_key.data() + num_atoms,
            m_solvent_atom_order.data());
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
   * was_overflown() check after the kernel runs (see below). Falls back
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
    // reported via was_overflown(), but the report went to io::messages
    // and nothing had displayed it yet, so it looked like silent data
    // corruption instead of the flagged error it actually was.
    const unsigned capacity = static_cast<unsigned>(std::ceil(per_block * num_blocks_a)) + num_blocks_a;
    return std::min(capacity, worst_case);
  }

  void check_candidate_overflow(gpu::TileVecT<gpu::Interaction_Tile> const & tiles,
                                 const char * which) {
    cudaDeviceSynchronize();
    if (tiles.was_overflown()) {
      io::messages.add(
        std::string("CUDA pairlist candidate build overflowed its ") + which +
        " capacity estimate -- candidates were dropped, the pairlist is "
        "wrong. This is a sizing-heuristic bug, not a real memory limit; "
        "increase estimate_candidate_capacity's margin.",
        "CUDA_Pairlist_Algorithm", io::message::error);
    }
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
    if (num_solute_blocks > 0) {
        const unsigned cap_ss = estimate_candidate_capacity(num_solute_blocks, num_solute_blocks, cutoff_d, vol, true);
        const unsigned cap_sv = num_solvent_blocks > 0
            ? estimate_candidate_capacity(num_solute_blocks, num_solvent_blocks, cutoff_d, vol, false)
            : 0;
        m_tiles.solute_candidates.reserve(cap_ss + cap_sv);
    }

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

    if (num_solute_blocks > 0) {
        check_candidate_overflow(m_tiles.solute_candidates, "solute_candidates");
    }

    if (num_solvent_blocks > 0) {
        // solvent - solvent (self-pairs, bi <= bj only)
        m_tiles.solvent_candidates.reserve(
            estimate_candidate_capacity(num_solvent_blocks, num_solvent_blocks, cutoff_d, vol, true));
        dim3 dimGrid((num_solvent_blocks + 15) / 16, (num_solvent_blocks + 15) / 16);
        gpu::find_block_candidates_kernel<B><<<dimGrid, dimBlock2D>>>(
            num_solvent_blocks, num_solvent_blocks,
            m_solvent_block_center.data(), m_solvent_block_radius.data(),
            m_solvent_block_center.data(), m_solvent_block_radius.data(),
            true, periodicity, cutoff, m_tiles.solvent_candidates.view());
        check_candidate_overflow(m_tiles.solvent_candidates, "solvent_candidates");
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
    if (m_tiles.solute_short.was_overflown()  || m_tiles.solute_long.was_overflown() ||
        m_tiles.solvent_short.was_overflown() || m_tiles.solvent_long.was_overflown()) {
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
                bool recompute_long)
{
    const unsigned num_atoms   = static_cast<unsigned>(m_force.size());
    const unsigned num_buckets = m_num_energy_groups * m_num_energy_groups;

    // IEEE-754 zero is the all-zero bit pattern -- cudaMemset is safe and
    // matches the existing convention (TileVecT::clear() zeroes its tile
    // array the same way).
    cudaMemset(m_force.data(), 0, sizeof(FPL3_TYPE) * num_atoms);
    cudaMemset(m_e_lj.data(),  0, sizeof(double) * num_buckets);
    cudaMemset(m_e_crf.data(), 0, sizeof(double) * num_buckets);

    // sync_pos_vel = false: prepare() already did this step's one real
    // resync; called from calculate_interactions(), always after
    // prepare().
    const math::CuVArray::View pos =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;
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
        m_force.data(), m_e_lj.data(), m_e_crf.data());

    gpu::launch_lj_crf_tiles(
        m_tiles.solvent_short.view(), m_solvent_atom_order.data(), num_solvent_atoms,
        nullptr, 0u,
        pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
        m_force.data(), m_e_lj.data(), m_e_crf.data());

    // Long-range: only recomputed on a classification/rebuild step
    // (recompute_long == true) -- otherwise m_longrange_force/
    // m_e_lj_long/m_e_crf_long are left untouched, holding whatever was
    // last computed here, matching nonbonded_set.cc's frozen
    // m_longrange_storage exactly. The zero must stay inside this
    // branch: zeroing unconditionally would wipe the frozen values a
    // non-rebuild step is supposed to reuse.
    if (recompute_long) {
        cudaMemset(m_longrange_force.data(), 0, sizeof(FPL3_TYPE) * num_atoms);
        cudaMemset(m_e_lj_long.data(),  0, sizeof(double) * num_buckets);
        cudaMemset(m_e_crf_long.data(), 0, sizeof(double) * num_buckets);

        gpu::launch_lj_crf_tiles(
            m_tiles.solute_long.view(), m_solute_atom_order.data(), num_solute_atoms,
            m_solvent_atom_order.data(), num_solvent_atoms,
            pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
            m_longrange_force.data(), m_e_lj_long.data(), m_e_crf_long.data());

        gpu::launch_lj_crf_tiles(
            m_tiles.solvent_long.view(), m_solvent_atom_order.data(), num_solvent_atoms,
            nullptr, 0u,
            pos, m_iac.data(), m_charge.data(), m_atom_energy_group.data(), lj, nb, boundary, box,
            m_longrange_force.data(), m_e_lj_long.data(), m_e_crf_long.data());
    }

    cudaDeviceSynchronize();

    // Forcefield::calculate_interactions() zeroes conf.current().force/
    // energies once before every Interaction in the sequence runs --
    // accumulate (+=), don't overwrite. Add the (possibly frozen-from-an-
    // earlier-call) long-range contribution in unconditionally, same as
    // nonbonded_set.cc's m_storage.force += m_longrange_storage.force.
    for (unsigned i = 0; i < num_atoms; ++i) {
        conf.current().force(i) += math::Vec(m_force[i].x, m_force[i].y, m_force[i].z);
        conf.current().force(i) += math::Vec(m_longrange_force[i].x,
                                              m_longrange_force[i].y,
                                              m_longrange_force[i].z);
    }

    // Direct per-[gi][gj] accumulation, same style as
    // nonbonded_innerloop.cc's CPU inner loop -- no scalar out-params.
    for (unsigned gi = 0; gi < m_num_energy_groups; ++gi) {
        for (unsigned gj = 0; gj < m_num_energy_groups; ++gj) {
            const unsigned k = gi * m_num_energy_groups + gj;
            conf.current().energies.lj_energy[gi][gj]  += m_e_lj[k]  + m_e_lj_long[k];
            conf.current().energies.crf_energy[gi][gj] += m_e_crf[k] + m_e_crf_long[k];
        }
    }
}
