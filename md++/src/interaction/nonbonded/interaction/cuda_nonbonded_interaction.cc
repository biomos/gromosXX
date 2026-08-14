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
 * GPU-native nonbonded interaction. See cuda_nonbonded_interaction.h.
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../math/boundary_checks.h"
#include "../../../math/gmath.h"

#include "../../../interaction/interaction.h"
#include "nonbonded_parameter.h"
#include "nonbonded_term.h"
#include "nonbonded_interaction.h"
#include "../pairlist/pairlist_algorithm.h"
#include "../pairlist/cuda_pairlist_algorithm.h"
#include "cuda_nonbonded_interaction.h"

#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

interaction::CUDA_Nonbonded_Interaction::CUDA_Nonbonded_Interaction(Pairlist_Algorithm * pa)
: Nonbonded_Interaction(pa)
{}

int interaction::CUDA_Nonbonded_Interaction::init(
                                      topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation & sim,
                                      std::ostream & os,
                                      bool quiet)
{
  if (!quiet)
    os << "CUDA NONBONDED INTERACTION\n";

  if (!math::boundary_check_cutoff(conf.current().box, conf.boundary_type,
            sim.param().pairlist.cutoff_long)) {
    io::messages.add("box is too small: not twice the cutoff!",
            "configuration", io::message::error);
    return 1;
  }

  // v1 scope (see this class's header comment): exactly one energy
  // group -- the tile kernel's energy reduction produces two flat
  // totals (e_lj, e_crf), not per-energy-group-pair buckets.
  if (topo.energy_groups().size() != 1) {
    io::messages.add(
      "CUDA_Nonbonded_Interaction only supports a single energy group so "
      "far (TILE_PAIRLIST_DESIGN.md §8/§9); the tile kernel's energy "
      "reduction produces one LJ and one CRF total, not a per-energy-"
      "group-pair breakdown.",
      "CUDA_Nonbonded_Interaction", io::message::error);
    return 1;
  }

  // v1 scope: no virial -- the tile kernel doesn't accumulate r (x) f.
  if (sim.param().pcouple.virial != math::no_virial) {
    io::messages.add(
      "CUDA_Nonbonded_Interaction does not compute a virial yet "
      "(TILE_PAIRLIST_DESIGN.md §8/§9); set PCOUPLE/VIRIAL to 0 "
      "(no virial) under accelerator = cuda.",
      "CUDA_Nonbonded_Interaction", io::message::error);
    return 1;
  }

  // v1 scope: no perturbation/EDS -- this class never builds a
  // Perturbed_Nonbonded_Set (or any Nonbonded_Set at all) and always
  // calls CUDA_Pairlist_Algorithm::update(), never update_perturbed().
  // Found the hard way: without this gate, a perturbed topology reaches
  // calculate_interactions() below and crashes (perturbed exclusions/
  // lambda-interpolated charges aren't what init()'s m_iac/m_charge
  // upload or the tile classification assume) instead of failing loudly
  // at init() the way every other unsupported-configuration gate in this
  // class does.
  if (sim.param().perturbation.perturbation || sim.param().eds.eds) {
    io::messages.add(
      "CUDA_Nonbonded_Interaction does not support perturbation or EDS "
      "yet (TILE_PAIRLIST_DESIGN.md §8/§9).",
      "CUDA_Nonbonded_Interaction", io::message::error);
    return 1;
  }

  // CUDA_Pairlist_Algorithm::init() has its own hard-error gate
  // (vacuum/rectangular boundary only, TILE_PAIRLIST_DESIGN.md §4.1) and
  // builds the per-atom iac/charge arrays this class's
  // compute_forces_energies() call needs every step.
  if (m_pairlist_algorithm->init(topo, conf, sim, os, quiet) != 0) {
    return 1;
  }

  // GPU-resident LJ parameter matrix, built once from m_parameter (the
  // CPU-side table create_nonbonded.cc already populated from the
  // topology file before calling this init()).
  m_gpu_lj.init(parameter());

  // Reaction-field constants, matching Nonbonded_Term::init's default
  // (eps=0, no coarse-graining) case exactly -- see lj_crf_tiles.h's doc
  // comment for the formulas these feed.
  interaction::Nonbonded_Term term;
  term.init(sim);
  m_nb.four_pi_eps_i   = static_cast<FPL_TYPE>(math::four_pi_eps_i);
  m_nb.crf_2cut3i      = static_cast<FPL_TYPE>(term.crf_2cut3i());
  m_nb.crf_cut         = static_cast<FPL_TYPE>(term.crf_cut(0));
  m_nb.cutoff_short_sq = static_cast<FPL_TYPE>(sim.param().pairlist.cutoff_short *
                                                sim.param().pairlist.cutoff_short);
  m_nb.cutoff_long_sq  = static_cast<FPL_TYPE>(sim.param().pairlist.cutoff_long *
                                                sim.param().pairlist.cutoff_long);

  m_initialized = true;

  if (!quiet)
    os << "END\n";
  return 0;
}

int interaction::CUDA_Nonbonded_Interaction::calculate_interactions(
                                      topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation & sim)
{
  DEBUG(7, "CUDA_Nonbonded_Interaction::calculate_interactions");

  if (!m_initialized) {
    // init() hard-errored (out-of-v1-scope configuration) -- see this
    // method's doc comment. No GPU state to touch safely; leave the
    // force/energy contribution at zero.
    return 1;
  }

  CUDA_Pairlist_Algorithm * pa =
      static_cast<CUDA_Pairlist_Algorithm *>(m_pairlist_algorithm);

  pa->prepare(topo, conf, sim); // every step: box-wrap + cog + GPU pos mirror refresh

  // Exact structural mirror of nonbonded_set.cc's pairlist_update check
  // (same operator, same operand order) -- real GROMOS twin-range: the
  // pairlist rebuild + classification, and the long-range force
  // recompute, only happen every skip_step steps; short-range forces
  // recompute every step regardless (see compute_forces_energies()).
  const bool pairlist_update = !(sim.steps() % sim.param().pairlist.skip_step);

  if (pairlist_update) {
    interaction::PairlistContainer dummy;
    dummy.resize(static_cast<unsigned>(topo.num_atoms()));
    pa->update(topo, conf, sim, dummy, 0, static_cast<unsigned>(topo.num_atoms()), 1);
  }

  double e_lj = 0.0, e_crf = 0.0;
  pa->compute_forces_energies(topo, conf, sim, m_gpu_lj.view(), m_nb,
                               pairlist_update, e_lj, e_crf);

  // Single energy group only (init()'s gate) -- everything lands in the
  // one (0, 0) bucket.
  conf.current().energies.lj_energy[0][0]  += e_lj;
  conf.current().energies.crf_energy[0][0] += e_crf;

  return 0;
}
