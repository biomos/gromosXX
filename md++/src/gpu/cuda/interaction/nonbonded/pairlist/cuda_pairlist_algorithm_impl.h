#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
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
       * put the chargegroups into the box. `topo_view`/`conf_view` are
       * fetched by prepare_cog() via sim.cuda().topology_view()/
       * configuration_view() (PLAN.md §3.2) -- passed in rather than
       * fetched here since this template has no simulation::Simulation&
       * of its own to call .cuda() on.
       */
      template<math::boundary_enum b>
      void _prepare_cog(configuration::Configuration & conf,
                        topology::Topology & topo,
                        gpu::Topology::View topo_view,
                        gpu::Configuration::View conf_view);

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
       * `conf_view` is fetched by reorder() via sim.cuda().
       * configuration_view() (PLAN.md §3.2) and passed in, same reason
       * as _prepare_cog() above.
       */
      template<math::boundary_enum b>
      void _atom_sort_key_atomic(configuration::Configuration & conf,
                                  unsigned num_atoms,
                                  gpu::Configuration::View conf_view);

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
       * TILE_PAIRLIST_DESIGN.md §10, part 2: decides whether the
       * (expensive) candidate rebuild is due this classification cycle,
       * decoupled from classify_tiles()'s own skip_step cadence via the
       * Verlet-buffer criterion: rebuild once
       * `2 * max_displacement_since_last_rebuild >= skin`, since two
       * atoms could in the worst case have closed the gap between them
       * by that much since the candidates were last built at
       * `cutoff_long + skin`.
       *
       * Always returns true if no candidate rebuild has ever happened
       * (m_candidates_built == false) or if skin == 0.0 -- the latter is
       * an explicit early-return, not left as an emergent property of
       * `2*0 >= 0`, so skin's degenerate case (rebuild every
       * classification cycle, today's exact pre-this-feature behavior)
       * is obviously correct on inspection rather than relying on a
       * floating-point comparison that could misbehave for a
       * near-zero-but-nonzero skin from a config round-trip.
       *
       * This check itself only runs at classification cadence (every
       * skip_step steps, not every step) since it's only ever called
       * from update() -- meaning true displacement is only ever known
       * "as of" the last classification check, up to skip_step steps
       * late. Accepted trade-off, not a silent gap: checking continuously
       * would defeat skip_step's whole performance purpose.
       */
      bool needs_candidate_rebuild(configuration::Configuration & conf,
                                    simulation::Simulation & sim);

      /**
       * Rebuilds the candidate tile lists (reorder() + build_candidates(),
       * unchanged), then snapshots the current GPU position mirror into
       * m_candidate_ref_pos so the next needs_candidate_rebuild() call
       * has something to diff against, and marks m_candidates_built.
       * This is what update() calls instead of reorder()+build_candidates()
       * directly, conditioned on needs_candidate_rebuild().
       */
      void rebuild_candidates(configuration::Configuration & conf,
                               topology::Topology & topo,
                               simulation::Simulation & sim);

      /**
       * Number of times rebuild_candidates() has actually run. Exposed
       * for the skin-buffer drift test (TILE_PAIRLIST_DESIGN.md §10 part
       * 2, PLAN.md §9.4) to confirm skin is actually reducing rebuild
       * frequency, not just numerically matching by coincidence (a bug
       * that made needs_candidate_rebuild() always return true would
       * still pass a pure force/energy comparison).
       */
      unsigned candidate_rebuild_count() const { return m_candidate_rebuild_count; }

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

      /**
       * The real production entry point (TILE_PAIRLIST_DESIGN.md §8/§9,
       * §10 twin-range cadence): runs gpu::lj_crf_tile_kernel over
       * solute_short/solvent_short every call, accumulating forces into a
       * device buffer and energies into two double accumulators. The
       * long-range contribution (solute_long/solvent_long) is only
       * recomputed when `recompute_long` is true (a classification/
       * rebuild step, per real GROMOS twin-range semantics -- see
       * nonbonded_set.cc's pairlist_update-gated longrange_storage) --
       * m_longrange_force/m_e_lj_long/m_e_crf_long hold whatever was last
       * computed and are added into every call's total regardless,
       * frozen in between. Adds (+=, not overwrites --
       * Forcefield::calculate_interactions zeroes conf.current().force/
       * energies once before every Interaction in the sequence runs) the
       * combined short+long result into conf.current().force and returns
       * the total LJ/CRF energies via e_lj/e_crf. Must be called after
       * classify_tiles(); m_iac/m_charge (built once in init(), topology
       * is static for a normal run) and m_force/m_e_lj/m_e_crf/
       * m_longrange_force/m_e_lj_long/m_e_crf_long (sized once in
       * init()) back this.
       */
      void compute_forces_energies(configuration::Configuration & conf,
                                    topology::Topology & topo,
                                    simulation::Simulation & sim,
                                    gpu::LJParamView lj,
                                    gpu::NbSimParams nb,
                                    bool recompute_long,
                                    double & e_lj,
                                    double & e_crf);

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

      /**
       * Per-atom integer atom code / charge, built once in init() (static
       * for a normal run, like the exclusion CSR in gpu::Topology) --
       * compute_forces_energies() reads these directly instead of going
       * through TopologyView every call.
       */
      gpu::cuvector<int> m_iac;
      gpu::cuvector<FPL_TYPE> m_charge;
      /**
       * Per-atom force accumulator and total LJ/CRF energy accumulators,
       * sized once in init(), re-zeroed at the top of every
       * compute_forces_energies() call.
       */
      gpu::cuvector<FPL3_TYPE> m_force;
      gpu::cuvector<double> m_e_lj;
      gpu::cuvector<double> m_e_crf;

      /**
       * Long-range (solute_long + solvent_long) per-atom force and total
       * energy, frozen between classification/rebuild steps -- mirrors
       * nonbonded_set.cc's m_longrange_storage exactly (real GROMOS
       * twin-range: the long-range force *value*, not just the pair
       * list, is held static between rebuilds, only recomputed on a
       * classification step). Sized once in init(); only ever written
       * inside compute_forces_energies() when its recompute_long
       * parameter is true -- never zeroed/touched otherwise, so a
       * non-rebuild step's call reuses whatever was last computed here.
       */
      gpu::cuvector<FPL3_TYPE> m_longrange_force;
      gpu::cuvector<double> m_e_lj_long;
      gpu::cuvector<double> m_e_crf_long;

      /**
       * Position snapshot at the time of the most recent CANDIDATE
       * rebuild (not the most recent classification) -- feeds
       * needs_candidate_rebuild()'s displacement check. Sized once in
       * init(); refreshed only inside rebuild_candidates().
       */
      math::CuVArray m_candidate_ref_pos;
      /**
       * Scratch for launch_max_displacement()'s per-block partial
       * maxima -- owned here so needs_candidate_rebuild() doesn't
       * allocate/free it on every classification-cadence call.
       */
      gpu::cuvector<FPL_TYPE> m_displacement_partial;
      /**
       * False until the first candidate rebuild ever happens -- there is
       * nothing to diff m_candidate_ref_pos against yet, so
       * needs_candidate_rebuild() must unconditionally return true.
       */
      bool m_candidates_built = false;
      unsigned m_candidate_rebuild_count = 0;
  };
}
