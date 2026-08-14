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
 * @file cuda_nonbonded_interaction.h
 * GPU-native nonbonded interaction (PLAN.md §10 step 9,
 * TILE_PAIRLIST_DESIGN.md §8/§9).
 *
 * Only ever included under USE_CUDA -- see create_nonbonded.cc, the only
 * call site (same convention as cuda_pairlist_algorithm.h).
 *
 * v1 scope, hard-errored in init() rather than silently producing wrong
 * numbers: exactly one energy group (Configuration::Energy::lj_energy/
 * crf_energy are per-energy-group-pair matrices; the tile kernel's
 * reduction only ever produces two flat totals, not per-group-pair
 * buckets -- generalizing that is real future work, not done here) and
 * no virial (sim.param().pcouple.virial must be math::no_virial; the tile
 * kernel doesn't accumulate r (x) f at all yet). Both restrictions inherit
 * from CUDA_Pairlist_Algorithm's own v1 scope (vacuum/rectangular
 * boundary only).
 *
 * Also does not implement the CPU twin-range performance optimization of
 * freezing long-range forces between pairlist rebuilds
 * (Nonbonded_Set::calculate_interactions holds solute_long/solvent_long
 * force contributions static between rebuilds, only solute_short/
 * solvent_short forces are recomputed every step) -- every call to
 * calculate_interactions() here fully re-runs prepare()+update()
 * (candidate build + classification) and recomputes forces/energies for
 * *all* four tile buckets from current positions, every step, regardless
 * of sim.param().pairlist.skip_step. This is numerically exact for any
 * single evaluation (right after a rebuild, "recompute fresh" and
 * "reuse what was frozen at this same rebuild" are the same numbers) and
 * therefore doesn't compromise a force/energy comparison test, but it
 * does not reproduce skip_step's performance characteristic, nor -- for
 * a multi-step trajectory where atoms drift across the cutoff_short/
 * cutoff_long boundary between what would have been rebuild steps -- the
 * exact per-step energies a real skip_step > 1 CPU run would produce.
 * Revisit once this needs to actually run production MD, not just be
 * numerically validated once.
 */

#pragma once

#include "nonbonded_interaction.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"

namespace interaction {

  class Pairlist_Algorithm;

  class CUDA_Nonbonded_Interaction : public Nonbonded_Interaction {
  public:
    explicit CUDA_Nonbonded_Interaction(Pairlist_Algorithm * pa);
    virtual ~CUDA_Nonbonded_Interaction() {}

    /**
     * Hard-errors (io::message::error, returns 1) for anything outside
     * v1 scope -- see this file's header comment -- rather than
     * delegating to Nonbonded_Interaction::init() (which would build CPU
     * Nonbonded_Set objects this class never uses). Still calls
     * m_pairlist_algorithm->init(...) (CUDA_Pairlist_Algorithm's own
     * boundary-scope gate) and builds the GPU-resident LJ parameter
     * matrix / reaction-field constants once.
     */
    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    /**
     * Real forces/energies, computed entirely via the tile pairlist +
     * gpu::lj_crf_tile_kernel (TILE_PAIRLIST_DESIGN.md §8) -- does not
     * use m_nonbonded_set (empty; never populated by this class's
     * init()) or any CPU Nonbonded_Innerloop/Outerloop code path.
     *
     * Returns 0 (no-op, force/energy contribution left at zero) without
     * touching any GPU state if init() didn't complete successfully --
     * found the hard way: some existing test harnesses (e.g.
     * aladip_cuda.t.cc) call Forcefield::init() without checking its
     * return value, so a hard-errored init() (out-of-v1-scope
     * configuration, e.g. perturbation) does not by itself stop
     * calculate_interactions() from being called afterward on
     * CUDA_Pairlist_Algorithm_Impl state (m_iac/m_charge/m_force/etc.)
     * that init() never allocated -- without this guard that's a crash
     * (CUDA calls on null/zero-capacity buffers), not just wrong numbers.
     */
    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

    /**
     * Not implemented for the GPU tile pipeline (single-pair evaluation
     * doesn't map onto "run the tile kernel"). Overridden solely to avoid
     * crashing: the base Nonbonded_Interaction::calculate_interaction()
     * unconditionally indexes m_nonbonded_set[0], which this class never
     * populates -- found via check_forcefield.cc's finite-difference
     * hessian check, which calls this on every Nonbonded_Interaction
     * subclass unconditionally, including one whose init() already
     * hard-errored (perturbation, out of v1 scope) and so was never
     * going to run this MD anyway.
     */
    virtual int calculate_interaction(topology::Topology & topo,
                                       configuration::Configuration & conf,
                                       simulation::Simulation & sim,
                                       unsigned int atom_i, unsigned int atom_j,
                                       math::Vec & force,
                                       double & e_lj, double & e_crf) {
      force = 0.0;
      e_lj = 0.0;
      e_crf = 0.0;
      return 1;
    }

  private:
    gpu::LJParams m_gpu_lj;
    gpu::NbSimParams m_nb;
    bool m_initialized = false;
  };

} // interaction
