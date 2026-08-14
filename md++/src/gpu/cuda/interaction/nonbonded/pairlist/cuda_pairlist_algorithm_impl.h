#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "interaction/nonbonded/pairlist/pairlist.h"

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
       * Sort solute and solvent ATOMS (separately) by their owning
       * chargegroup's Morton cell key, chunk each sorted sequence into
       * fixed 32-wide blocks, and compute each block's bounding sphere
       * from the 32 atoms' own positions. TILE_PAIRLIST_DESIGN.md §3
       * steps 1-3. Atom-indexed, not chargegroup-indexed -- see that
       * section's "second correction" for why (exclusions are atom-pair
       * granularity, chargegroup-cutoff's cog-cog test is a per-pair
       * lookup done later in classify_tiles(), not a property of the
       * blocks themselves).
       */
      void reorder(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      /**
       * atomic_cutoff mode's atom sort key: computed directly from each
       * atom's own position (a locally box-wrapped copy for the cell
       * lookup only, TILE_PAIRLIST_DESIGN.md §4.2/step 7), not inherited
       * from a chargegroup. Called from reorder() instead of
       * atom_sort_key_kernel when sim.param().pairlist.atomic_cutoff.
       */
      template<math::boundary_enum b>
      void _atom_sort_key_atomic(configuration::Configuration & conf,
                                  unsigned num_atoms);

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
       * Exclusion + short/long classification (TILE_PAIRLIST_DESIGN.md §3
       * step 5, §4.2/step 7): resolves every candidate tile's atom pairs,
       * applies exclusions, and buckets survivors into m_tiles.solute_short/
       * solute_long/solvent_short/solvent_long, using either chargegroup-
       * cog-cog distance or atom-atom distance depending on
       * sim.param().pairlist.atomic_cutoff. Must be called after
       * build_candidates().
       */
      void classify_tiles(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      template<math::boundary_enum b>
      void _classify_tiles(configuration::Configuration & conf,
                           topology::Topology & topo,
                           simulation::Simulation & sim);

      /**
       * Tiles from the most recent build_candidates()/classify_tiles()
       * call. Exposed for the pairlist-equivalence test (TILE_PAIRLIST_
       * DESIGN.md §5), which doesn't exist yet, and eventually a force
       * kernel (out of scope here).
       */
      const gpu::TileContainer & tiles() const { return m_tiles; }

      /**
       * Convert the classified tiles into a CPU-comparable
       * interaction::PairlistContainer, unpacking each tile's mask bits
       * back into (atom_i, atom_j) pairs via the same row/col atom-order
       * arrays used to build them. For the pairlist-equivalence test
       * (TILE_PAIRLIST_DESIGN.md §5/§6) only -- nothing in the production
       * pipeline needs the CPU-shaped format, and this is not a
       * performant way to get there (it's a host-side, per-bit unpack).
       */
      interaction::PairlistContainer to_pairlist_container(topology::Topology & topo) const;

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
       * Per-chargegroup representative (box-wrapped) position: true cog
       * for solute chargegroups, first-atom position for solvent (see
       * Periodicity::prepare_chargegroup) -- sized to *all* chargegroups.
       * Used by classify_tiles() for the chargegroup-cutoff cog-cog
       * distance test; no longer used for block bounding spheres (those
       * are atom-indexed now, see m_solute_atom_order below).
       */
      math::CuVArray m_cg_cog;
      /**
       * chargegroup cell indices array.
       */
      gpu::cuvector<ushort4> m_cg_cells;
      /**
       * flat copy of m_cg_cells[i].w (Morton cell index), one entry per
       * chargegroup -- feeds atom_sort_key_kernel (every atom inherits
       * its owning chargegroup's key) and is itself sorted in place by
       * reorder() to indirectly not need a separate atom-level Morton
       * computation.
       */
      gpu::cuvector<unsigned> m_cg_sort_key;
      /**
       * flat copy of m_cg_sort_key indexed by ATOM (atom_sort_key[a] =
       * cg_sort_key[owning_chargegroup(a)]), one entry per atom. This is
       * what actually gets thrust::sort_by_key'd in reorder() -- blocks
       * are atom-indexed (TILE_PAIRLIST_DESIGN.md §3's "second
       * correction"), so the sort key needs atom granularity too.
       */
      gpu::cuvector<unsigned> m_atom_sort_key;

      /**
       * Block-sorted permutations: solute/solvent block-sorted position ->
       * global atom index (solvent values are offset by
       * num_solute_atoms). Built by reorder().
       */
      gpu::cuvector<unsigned> m_solute_atom_order;
      gpu::cuvector<unsigned> m_solvent_atom_order;

      /**
       * Per-block bounding sphere (center + radius), one entry per
       * ceil(count / gpu::BLOCK_SIZE) block, computed directly from the
       * 32 atoms' own positions. Built by reorder(), consumed by
       * build_candidates().
       */
      gpu::cuvector<FPL3_TYPE> m_solute_block_center;
      gpu::cuvector<FPL_TYPE>  m_solute_block_radius;
      gpu::cuvector<FPL3_TYPE> m_solvent_block_center;
      gpu::cuvector<FPL_TYPE>  m_solvent_block_radius;

      /**
       * Candidate and short/long-classified tiles. See tiles() accessor.
       */
      gpu::TileContainer m_tiles;
  };
}
