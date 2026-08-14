/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 *
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file cuda_pairlist_algorithm.h
 * CUDA accelerated pairlist algorithm
 *
 * Only ever included under USE_CUDA -- see create_nonbonded.cc, the only
 * call site. Per PLAN.md D5 there is exactly one GPU-native pairlist and
 * it is never dispatched through make_algorithm's Backend-template
 * mechanism, so this class deliberately isn't a Backend template either
 * (it used to be; that bought no genericity since the <gpuBackend>
 * specialization shared nothing with the primary/cpuBackend template
 * body -- dropped in favor of this plain class + the standard .h/.cc/.cu
 * split already used for every other GPU-only feature).
 */

#pragma once

#include "gpu/cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"

namespace math
{
  template<math::boundary_enum>
  class Periodicity;
}

namespace interaction
{
  /**
   * @class CUDA_Pairlist_Algorithm
   * create an atomic pairlist on GPU
   * with a chargegroup based or atom
   * based cut-off criterion.
   */
  class CUDA_Pairlist_Algorithm : public Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    CUDA_Pairlist_Algorithm();

    /**
     * destructor.
     */
    virtual ~CUDA_Pairlist_Algorithm() {}

    /**
     * init. Hard-errors if the boundary type isn't one of the ones the
     * block/candidate build actually supports (TILE_PAIRLIST_DESIGN.md
     * §4.1: vacuum + rectangular only for v1 -- gpu::Periodicity::
     * set_cell_size has no real triclinic/truncoct implementation yet).
     * Both chargegroup-cutoff and atomic-cutoff (sim.param().pairlist.
     * atomic_cutoff) are supported (§4.2/step 7); no boundary restriction
     * beyond the one above applies to either.
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * prepare the pairlist(s).
     */
    virtual int prepare(topology::Topology & topo,
                        configuration::Configuration & conf,
                        simulation::Simulation &sim);

    /**
     * update the pairlist(s).
     */
    virtual void update(topology::Topology & topo,
                        configuration::Configuration & conf,
                        simulation::Simulation &sim,
                        interaction::PairlistContainer &pairlist,
                        unsigned int begin, unsigned int end,
                        unsigned int stride);

    virtual void update_perturbed(topology::Topology & topo,
                                  configuration::Configuration & conf,
                                  simulation::Simulation & sim,
                                  interaction::PairlistContainer & pairlist,
                                  interaction::PairlistContainer & perturbed_pairlist,
                                  unsigned int begin, unsigned int end,
                                  unsigned int stride) {
      io::messages.add(
        "update_perturbed is not supported by CUDA_Pairlist_Algorithm",
        "CUDA_Pairlist_Algorithm", io::message::error);
    };

    /**
     * For the pairlist-equivalence test (TILE_PAIRLIST_DESIGN.md §5/§6)
     * only -- converts the tiles built by the most recent update() call
     * into a CPU-comparable interaction::PairlistContainer. Not part of
     * the production pipeline; update()'s own `pairlist` output
     * parameter stays the intentional dummy (PAIRLIST_PLAN.md §5(A)).
     */
    interaction::PairlistContainer to_pairlist_container(topology::Topology & topo) const {
      return m_impl.to_pairlist_container(topo);
    }

    /**
     * The real production entry point (TILE_PAIRLIST_DESIGN.md §8/§9,
     * PLAN.md §10 step 9/twin-range cadence): runs the LJ + reaction-
     * field tile kernel over the tiles built by the most recent
     * update() call and accumulates the result into conf.current().
     * force (+=) and conf.current().energies.lj_energy/crf_energy's
     * per-energy-group-pair matrix (+=), same accumulation style as
     * nonbonded_innerloop.cc's CPU inner loop. Called by
     * CUDA_Nonbonded_Interaction, never directly by anything
     * CPU-pairlist-shaped -- unlike update()'s own `pairlist`
     * parameter, this is where the real numbers come from.
     *
     * `recompute_long`: real GROMOS twin-range only recomputes the
     * long-range (solute_long/solvent_long) contribution on a
     * classification/rebuild step -- pass true only when the caller
     * also called update() this same cycle (see CUDA_Nonbonded_
     * Interaction::calculate_interactions()'s pairlist_update check).
     * Short-range is always recomputed regardless.
     */
    void compute_forces_energies(topology::Topology & topo,
                                  configuration::Configuration & conf,
                                  simulation::Simulation & sim,
                                  gpu::LJParamView lj,
                                  gpu::NbSimParams nb,
                                  bool recompute_long) {
      m_impl.compute_forces_energies(conf, topo, sim, lj, nb, recompute_long);
    }

    /**
     * Test-support only (TILE_PAIRLIST_DESIGN.md §10 part 2's drift
     * test) -- see CUDA_Pairlist_Algorithm_Impl::candidate_rebuild_count().
     */
    unsigned candidate_rebuild_count() const { return m_impl.candidate_rebuild_count(); }

  private:
    CUDA_Pairlist_Algorithm_Impl m_impl;
  };
} // interaction
