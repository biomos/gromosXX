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
 * @file cuda_nonbonded_interaction.cc
 * CUDA_Nonbonded_Interaction: vanilla GPU LJ + CRF nonbonded forces.
 *
 * Pair list construction is delegated to the CPU Standard_Pairlist_Algorithm
 * which handles all exclusions, chargegroups, and periodic images correctly.
 * The resulting flat pair list is then processed by GPU force kernels.
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"
#include "simulation/parameter.h"

#include "interaction/interaction.h"
#include "interaction/interaction_types.h"
#include "interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "interaction/nonbonded/interaction/cuda_nonbonded_interaction.h"

#include "gpu/cuda/interaction/nonbonded/kernels/nb_kernels.h"
#include "gpu/cuda/utils.h"

#include "math/periodicity.h"
#include "math/boundary_checks.h"
#include "util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

// ─────────────────────────────────────────────────────────────────────────────

interaction::CUDA_Nonbonded_Interaction::CUDA_Nonbonded_Interaction(
    CUDA_Pairlist_Algorithm<util::gpuBackend>* pa)
    : Nonbonded_Interaction(pa)
{
}

interaction::CUDA_Nonbonded_Interaction::~CUDA_Nonbonded_Interaction()
{
    delete m_pairlist_algorithm;
}

// ─────────────────────────────────────────────────────────────────────────────
// build_flat_pairlist: convert CPU PairlistContainer → flat GPU cuvector
// ─────────────────────────────────────────────────────────────────────────────

void interaction::CUDA_Nonbonded_Interaction::build_flat_pairlist(
    const PairlistContainer& pl)
{
    m_gpu_pairs.clear();

    auto append = [&](const Pairlist& src) {
        for (unsigned i = 0; i < static_cast<unsigned>(src.size()); ++i) {
            for (unsigned j : src[i]) {
                m_gpu_pairs.push_back(make_uint2(i, j));
            }
        }
    };

    append(pl.solute_short);
    append(pl.solute_long);
    append(pl.solvent_short);
    append(pl.solvent_long);
}

// ─────────────────────────────────────────────────────────────────────────────
// init
// ─────────────────────────────────────────────────────────────────────────────

int interaction::CUDA_Nonbonded_Interaction::init(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim,
    std::ostream&                 os,
    bool                          quiet)
{
    if (!quiet) os << "NONBONDED INTERACTION (CUDA)\n";

    if (!math::boundary_check_cutoff(conf.current().box, conf.boundary_type,
                                     sim.param().pairlist.cutoff_long)) {
        io::messages.add("box too small: not twice the long cutoff",
                         "CUDA_Nonbonded_Interaction", io::message::error);
        return 1;
    }

    // ── Initialise CPU pairlist algorithm ────────────────────────────────────
    if (m_std_pairlist_alg.init(topo, conf, sim, os, quiet))
        return 1;

    m_pairlist.resize(topo.num_atoms());

    // ── Build GPU LJ parameter matrix ────────────────────────────────────────
    m_lj_params.init(m_parameter);

    // ── Precompute reaction-field constants (mirrors Nonbonded_Term::init) ───
    {
        const auto& p   = sim.param().nonbonded;
        const double rc = p.rf_cutoff;

        double crf = 0.0;
        if (rc > 0.0) {
            if (p.rf_epsilon == 0.0) {
                crf = -1.0;
            } else {
                const double kr  = p.rf_kappa * rc;
                const double kr2 = kr * kr;
                crf = (2.0 * (p.epsilon - p.rf_epsilon) * (1.0 + kr)
                       - p.rf_epsilon * kr2)
                    / ((p.epsilon + 2.0 * p.rf_epsilon) * (1.0 + kr)
                       + p.rf_epsilon * kr2);
            }
            const double cut3i = 1.0 / (rc * rc * rc);
            m_nb_params.crf_2cut3i = static_cast<FPL_TYPE>(crf * cut3i / 2.0);
            m_nb_params.crf_cut    = static_cast<FPL_TYPE>((1.0 - crf / 2.0) / rc);
        }

        m_nb_params.four_pi_eps_i =
            static_cast<FPL_TYPE>(math::four_pi_eps_i);
        m_nb_params.cutoff_short_sq =
            static_cast<FPL_TYPE>(sim.param().pairlist.cutoff_short *
                                  sim.param().pairlist.cutoff_short);
        m_nb_params.cutoff_long_sq =
            static_cast<FPL_TYPE>(sim.param().pairlist.cutoff_long *
                                  sim.param().pairlist.cutoff_long);
    }

    // ── Reserve GPU pair list (rough estimate) ────────────────────────────────
    m_gpu_pairs.reserve(topo.num_atoms() * 80u);

    // ── Initial full data upload to GPU ──────────────────────────────────────
    conf.copy_to_gpu();
    // Ensure GPU force array is the right size
    conf.get_gpu_raw_ptrs(); // triggers lazy GPU Configuration resize if needed

    if (!quiet) {
        os << "\tGPU force kernel  : LJ + reaction-field CRF\n";
        os << "\tLJ types          : " << m_lj_params.num_types << "\n";
        os << "\tcrf_2cut3i        : " << m_nb_params.crf_2cut3i << "\n";
        os << "\tcrf_cut           : " << m_nb_params.crf_cut    << "\n";
        os << "END\n";
    }

    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// calculate_interactions
// ─────────────────────────────────────────────────────────────────────────────

int interaction::CUDA_Nonbonded_Interaction::calculate_interactions(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    DEBUG(4, "CUDA_Nonbonded_Interaction::calculate_interactions");
    m_timer.start(sim);

    int steps = sim.param().multistep.steps;
    if (steps == 0) steps = 1;

    if ((sim.steps() % steps) == 0) {

        // ── Step 1: CPU pairlist construction ────────────────────────────────
        const bool pairlist_update =
            !(sim.steps() % sim.param().pairlist.skip_step);

        if (pairlist_update) {
            if (m_std_pairlist_alg.prepare(topo, conf, sim))
                return 1;
            m_std_pairlist_alg.update(topo, conf, sim,
                                      m_pairlist, 0,
                                      topo.num_atoms(), 1);
        }

        // ── Step 2: Convert pairlist to flat GPU format ───────────────────────
        build_flat_pairlist(m_pairlist);

        // ── Step 3: Upload current positions to GPU ───────────────────────────
        conf.copy_pos_vel_to_gpu();

        // ── Step 4: Zero GPU force array ──────────────────────────────────────
        auto ptrs = conf.get_gpu_raw_ptrs();
        const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
        gpu::launch_zero_forces(ptrs.force, num_atoms);

        // ── Step 5: LJ-CRF force kernel ───────────────────────────────────────
        gpu::launch_lj_crf_forces(
            ptrs.pos,
            ptrs.force,
            topo.get_gpu_view().iac,
            topo.get_gpu_view().charge,
            m_gpu_pairs.data(),
            static_cast<unsigned>(m_gpu_pairs.size()),
            m_lj_params.view(),
            m_nb_params,
            conf.boundary_type,
            conf.current().box);

        // ── Step 6: Sync GPU forces back to CPU ───────────────────────────────
        conf.copy_forces_from_gpu();
    }

    DEBUG(6, "CUDA_Nonbonded_Interaction::calculate_interactions done");
    m_timer.stop();
    return 0;
}
