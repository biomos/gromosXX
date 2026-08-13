#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/pairlist/tile.h"

namespace interaction {
  /**
   * @class CUDA_Pairlist_Algorithm_Impl
   * GPU-native implementation backing CUDA_Pairlist_Algorithm.
   *
   * Not a Backend-templated pImpl: CUDA_Pairlist_Algorithm is only ever
   * compiled/included under USE_CUDA (its single call site,
   * create_nonbonded.cc, guards the #include the same way), so there is
   * no CPU counterpart to keep in sync with. Per PLAN.md D5, there is
   * exactly one GPU-native pairlist and it is never dispatched through
   * make_algorithm's Backend-template mechanism -- the earlier
   * Backend-templated shape bought no genericity (the <gpuBackend>
   * specialization shared nothing with the primary template body) and
   * was dropped for this plain class + the standard .h/.cc/.cu split.
   */
  class CUDA_Pairlist_Algorithm_Impl {
    public:
      CUDA_Pairlist_Algorithm_Impl();

      int init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os = std::cout,
        bool quiet = false);

      void set_cutoff(double const cutoff_short, double const cutoff_long);

      void prepare_cog(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      /**
       * put the chargegroups into the box
       */
      template<math::boundary_enum b>
      void _prepare_cog(configuration::Configuration & conf,
                        topology::Topology & topo);

      /**
       * Sort solute and solvent chargegroups (separately) by their
       * Morton cell key, chunk each sorted sequence into fixed 32-wide
       * blocks, and compute each block's bounding sphere.
       * TILE_PAIRLIST_DESIGN.md §3 steps 2-3.
       */
      void reorder(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      /**
       * Block-pair candidate search (TILE_PAIRLIST_DESIGN.md §3 step 4):
       * test every solute-solute, solute-solvent, and solvent-solvent
       * block-pair's bounding-sphere distance against cutoff_long + skin,
       * and push survivors into m_tiles.solute_candidates /
       * m_tiles.solvent_candidates. Must be called after reorder().
       */
      void build_candidates(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      template<math::boundary_enum b>
      void _build_candidates(configuration::Configuration & conf,
                              topology::Topology & topo,
                              simulation::Simulation & sim);

      /**
       * Candidate tiles from the most recent build_candidates() call.
       * Not yet consumed by anything -- exposed for the classification
       * pass (TILE_PAIRLIST_DESIGN.md §3 step 5) and the pairlist-
       * equivalence test (§5), neither of which exist yet.
       */
      const gpu::TileContainer & tiles() const { return m_tiles; }

    protected:
      /**
       * squared shortrange cutoff.
       */
      double m_cutoff_short_2;
      /**
       * squared longrange cutoff.
       */
      double m_cutoff_long_2;
      /**
       * longrange cutoff.
       */
      double m_cutoff_long;
      /**
       * shortrange cutoff.
       */
      double m_cutoff_short;

    private:
      /**
       * chargegroup center of geometry array.
       */
      math::CuVArray m_cg_cog;
      /**
       * chargegroup cell indices array.
       */
      gpu::cuvector<ushort4> m_cg_cells;
      /**
       * flat copy of m_cg_cells[i].w (Morton cell index), one entry per
       * chargegroup -- feeds a Thrust sort-by-key in the block-build step
       * (TILE_PAIRLIST_DESIGN.md §3 step 2), which needs a plain array,
       * not one field of a ushort4.
       */
      gpu::cuvector<unsigned> m_cg_sort_key;

      /**
       * Block-sorted permutations: solute/solvent block-sorted position ->
       * global chargegroup index (solvent values are offset by
       * num_solute_chargegroups, since m_cg_cog/m_cg_cells/m_cg_sort_key
       * are indexed globally). Built by reorder().
       */
      gpu::cuvector<unsigned> m_solute_order;
      gpu::cuvector<unsigned> m_solvent_order;

      /**
       * Per-block bounding sphere (center + radius), one entry per
       * ceil(count / gpu::BLOCK_SIZE) block. Built by reorder(), consumed
       * by build_candidates().
       */
      gpu::cuvector<FPL3_TYPE> m_solute_block_center;
      gpu::cuvector<FPL_TYPE>  m_solute_block_radius;
      gpu::cuvector<FPL3_TYPE> m_solvent_block_center;
      gpu::cuvector<FPL_TYPE>  m_solvent_block_radius;

      /**
       * Candidate (and eventually short/long-classified) tiles. See
       * tiles() accessor.
       */
      gpu::TileContainer m_tiles;
  };
}
