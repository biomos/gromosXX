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

    if (!sim.param().pairlist.atomic_cutoff){
        conf.copy_to_gpu();
        const size_t num_solute_cg = topo.num_solute_chargegroups();
        const size_t num_cg = topo.num_chargegroups();
        m_cg_cog.resize(num_solute_cg);
        m_cg_cells.resize(num_cg);
        m_cg_sort_key.resize(num_cg);
        SPLIT_BOUNDARY(_prepare_cog, conf, topo);
    }
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_prepare_cog<B>(
                                    configuration::Configuration & conf,
                                    topology::Topology & topo) {
    DEBUG(10, "putting chargegroups into box");

    const unsigned num_cg = static_cast<unsigned>(topo.num_chargegroups());
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((num_cg + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();
    gpu::Periodicity<B> periodicity(conf.current().box);
    periodicity.set_cell_size(m_cutoff_long);
    gpu::prepare_cog_kernel<<<dimGrid, dimBlock>>>(topo.get_gpu_view(),
                                                    conf.get_gpu_view(),
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

    const gpu::Topology::View topo_view = topo.get_gpu_view();

    // Every atom inherits its owning chargegroup's Morton cell key
    // (TILE_PAIRLIST_DESIGN.md §3 step 1) -- keeps chargegroup members
    // spatially adjacent after sorting, without a separate per-atom
    // Morton computation.
    m_atom_sort_key.resize(num_atoms);
    {
        dim3 dimBlock(NUM_THREADS_PER_BLOCK);
        dim3 dimGrid((num_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
        gpu::atom_sort_key_kernel<<<dimGrid, dimBlock>>>(
            topo_view.chargegroup, num_cg, num_atoms,
            m_cg_sort_key.data(), m_atom_sort_key.data());
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

    // Bounding spheres come directly from each atom's own (box-wrapped)
    // position now, not a chargegroup cog -- m_cg_cog is reserved for
    // classify_tiles()'s chargegroup-cutoff distance test.
    const math::CuVArray::View pos = conf.get_gpu_view().current().pos;

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
    const unsigned capacity = static_cast<unsigned>(per_block * num_blocks_a) + num_blocks_a;
    return std::min(capacity, worst_case);
  }

  void check_candidate_overflow(gpu::TileVecT<gpu::Interaction_Tile> & tiles,
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
            true, periodicity, cutoff, m_tiles.solute_candidates);
    }

    if (num_solute_blocks > 0 && num_solvent_blocks > 0) {
        // solute - solvent (every pair; column side is solvent-order)
        dim3 dimGrid((num_solute_blocks + 15) / 16, (num_solvent_blocks + 15) / 16);
        gpu::find_block_candidates_kernel<B><<<dimGrid, dimBlock2D>>>(
            num_solute_blocks, num_solvent_blocks,
            m_solute_block_center.data(), m_solute_block_radius.data(),
            m_solvent_block_center.data(), m_solvent_block_radius.data(),
            false, periodicity, cutoff, m_tiles.solute_candidates);
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
            true, periodicity, cutoff, m_tiles.solvent_candidates);
        check_candidate_overflow(m_tiles.solvent_candidates, "solvent_candidates");
    }
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
    const gpu::Topology::View topo_view = topo.get_gpu_view();

    dim3 dimBlock2D(gpu::BLOCK_SIZE, gpu::BLOCK_SIZE);

    if (num_solute_candidates > 0) {
        gpu::classify_tiles_kernel<B><<<num_solute_candidates, dimBlock2D>>>(
            m_tiles.solute_candidates,
            m_solute_atom_order.data(), static_cast<unsigned>(m_solute_atom_order.size()),
            m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
            m_cg_cog.view(), topo_view, periodicity,
            cutoff_short2, cutoff_long2,
            m_tiles.solute_short, m_tiles.solute_long);
    }

    if (num_solvent_candidates > 0) {
        gpu::classify_tiles_kernel<B><<<num_solvent_candidates, dimBlock2D>>>(
            m_tiles.solvent_candidates,
            m_solvent_atom_order.data(), static_cast<unsigned>(m_solvent_atom_order.size()),
            nullptr, 0u,
            m_cg_cog.view(), topo_view, periodicity,
            cutoff_short2, cutoff_long2,
            m_tiles.solvent_short, m_tiles.solvent_long);
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
