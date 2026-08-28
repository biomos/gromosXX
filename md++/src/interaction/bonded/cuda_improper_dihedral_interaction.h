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
 * @file cuda_improper_dihedral_interaction.h
 * GPU-native improper dihedral interaction (PLAN.md §10 step 14,
 * bonded forces). Only ever included under USE_CUDA -- see
 * create_bonded.cc, the only call site. Same shape as
 * cuda_quartic_bond_interaction.h/cuda_angle_interaction.h: static term
 * list uploaded once, one-thread-per-term kernel every call, no
 * pairlist.
 *
 * v1 scope, hard-errored in init(): no perturbation, no GAMD (matching
 * every other CUDA gate in this codebase).
 */

#pragma once

#include "../interaction.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {

  class CUDA_Improper_Dihedral_Interaction : public Interaction {
  public:
    CUDA_Improper_Dihedral_Interaction() : Interaction("ImproperDihedral") {}
    virtual ~CUDA_Improper_Dihedral_Interaction();

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

    // See cuda_angle_interaction.h's identical override for why.
    virtual bool needs_fresh_cpu_force() const override { return false; }
    virtual bool is_gpu_native() const override { return true; }

  private:
    gpu::cuvector<unsigned> m_dihedral_i;
    gpu::cuvector<unsigned> m_dihedral_j;
    gpu::cuvector<unsigned> m_dihedral_k;
    gpu::cuvector<unsigned> m_dihedral_l;
    gpu::cuvector<unsigned> m_dihedral_type;
    gpu::cuvector<unsigned> m_atom_energy_group;
    gpu::cuvector<FPL_TYPE> m_K;
    gpu::cuvector<FPL_TYPE> m_q0;
    gpu::cuvector<double>   m_improper_energy;
    gpu::cuvector<double>   m_virial;
    unsigned m_num_dihedrals = 0;
    bool m_initialized = false;
    bool m_energy_registered = false;

    // Own stream: force/virial/energy are all written directly into
    // the GPU-resident mirror (atomicAdd) with no host sync at all --
    // m_improper_energy/m_virial are private per-call scratch for the
    // kernel to write into before that on-device merge, never touched
    // from the host.
    cudaStream_t m_stream = 0;
  };

} // namespace interaction
