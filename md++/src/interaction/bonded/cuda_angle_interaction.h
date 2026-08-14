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
 * @file cuda_angle_interaction.h
 * GPU-native harmonic (cosine) bond-angle interaction (PLAN.md §10 step
 * 14, bonded forces). Only ever included under USE_CUDA -- see
 * create_bonded.cc, the only call site. Same shape as
 * cuda_quartic_bond_interaction.h: static term list uploaded once,
 * one-thread-per-term kernel every call, no pairlist.
 *
 * v1 scope, hard-errored in init(): no perturbation, no GAMD (matching
 * every other CUDA gate in this codebase).
 */

#pragma once

#include "../interaction.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {

  class CUDA_Angle_Interaction : public Interaction {
  public:
    CUDA_Angle_Interaction() : Interaction("Angle") {}
    virtual ~CUDA_Angle_Interaction() {}

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

  private:
    gpu::cuvector<unsigned> m_angle_i;
    gpu::cuvector<unsigned> m_angle_j;
    gpu::cuvector<unsigned> m_angle_k;
    gpu::cuvector<unsigned> m_angle_type;
    gpu::cuvector<unsigned> m_atom_energy_group;
    gpu::cuvector<FPL_TYPE> m_K;
    gpu::cuvector<FPL_TYPE> m_cos0;
    gpu::cuvector<double>   m_angle_energy;
    gpu::cuvector<double>   m_virial;
    unsigned m_num_angles = 0;
    bool m_initialized = false;
  };

} // namespace interaction
