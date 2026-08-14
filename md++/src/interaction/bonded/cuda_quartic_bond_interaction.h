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
 * @file cuda_quartic_bond_interaction.h
 * GPU-native quartic bond interaction (PLAN.md §10 step 14, bonded
 * forces). Only ever included under USE_CUDA -- see create_bonded.cc,
 * the only call site (same convention as cuda_nonbonded_interaction.h).
 *
 * Unlike nonbonded, bonded terms have no pairlist/spatial search --
 * the term list (atom index pairs + type index) is static, read once
 * from topology, and small (aladip: ~80 bonds), so this class uploads
 * it once in init() and does a direct one-thread-per-term kernel launch
 * every calculate_interactions() call (no shared-memory bucketing,
 * unlike the nonbonded tile kernel's per-block-of-many-pairs shape).
 *
 * v1 scope, hard-errored in init(): no perturbation, no GAMD (matching
 * every other CUDA gate in this codebase).
 */

#pragma once

#include "../interaction.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {

  class CUDA_Quartic_Bond_Interaction : public Interaction {
  public:
    CUDA_Quartic_Bond_Interaction() : Interaction("QuarticBond") {}
    virtual ~CUDA_Quartic_Bond_Interaction() {}

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

  private:
    gpu::cuvector<unsigned> m_bond_i;
    gpu::cuvector<unsigned> m_bond_j;
    gpu::cuvector<unsigned> m_bond_type;
    gpu::cuvector<unsigned> m_atom_energy_group;
    gpu::cuvector<FPL_TYPE> m_K;
    gpu::cuvector<FPL_TYPE> m_r0;
    gpu::cuvector<double>   m_bond_energy;
    gpu::cuvector<double>   m_virial;
    unsigned m_num_bonds = 0;
    bool m_initialized = false;
  };

} // namespace interaction
