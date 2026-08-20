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
 * @file cuda_position_restraint_interaction.h
 * GPU-native position-restraint interaction. Only ever included under
 * USE_CUDA -- see create_special.cc, the only call site (same
 * convention as cuda_quartic_bond_interaction.h/cuda_angle_interaction.h).
 *
 * Unlike the bonded terms, the restraint list is a *sparse* subset of
 * atoms (topo.position_restraints(), not every atom), so the term list
 * (atom index + reference position + bfactor scale) is compact -- sized
 * to the number of restraints, not the number of atoms -- and, like the
 * bonded terms, uploaded once in init() (static for a normal run) with
 * one thread per restraint term every calculate_interactions() call.
 *
 * v1 scope: only NTPOR 1/2 (posrest_on/posrest_bfactor, soft restraining)
 * reach this class at all -- create_special.cc only builds a Position_
 * Restraint_Interaction (CPU or GPU) for those two parameter values.
 * NTPOR=3 (posrest_const, position *constraining*) is a completely
 * different mechanism (algorithm::Position_Constraints), not this class.
 * No perturbation gate: unlike every other CUDA_* class in this codebase,
 * Position_Restraint_Interaction has no perturbed counterpart anywhere
 * and its CPU calculate_interactions() never branches on
 * sim.param().perturbation.perturbation, so there is nothing to gate.
 */

#pragma once

#include "../interaction.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {

  class CUDA_Position_Restraint_Interaction : public Interaction {
  public:
    CUDA_Position_Restraint_Interaction() : Interaction("PositionRestraint") {}
    virtual ~CUDA_Position_Restraint_Interaction();

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

    // See cuda_angle_interaction.h's identical overrides for why: writes
    // force directly into the GPU-resident mirror via mark_gpu_dirty(
    // MIRROR_FORCE), must NOT trigger Forcefield::calculate_interactions
    // ()'s default pre-call flush.
    virtual bool needs_fresh_cpu_force() const override { return false; }
    virtual bool is_gpu_native() const override { return true; }

  private:
    gpu::cuvector<unsigned> m_seq;
    gpu::cuvector<FPL3_TYPE> m_ref;
    // 1/bf per restraint term (bf = bfactor(seq) under posrest_bfactor,
    // 1.0 otherwise) -- resolved once here so the kernel never branches
    // on posrest mode; see position_restraint_kernels.h's doc comment.
    gpu::cuvector<FPL_TYPE> m_inv_bf_scale;
    gpu::cuvector<unsigned> m_atom_energy_group;
    gpu::cuvector<double>   m_posrest_energy;
    unsigned m_num_restraints = 0;
    bool m_initialized = false;

    // Own stream: force is written directly into the GPU-resident
    // mirror (no sync at all); energy stays on a small private buffer
    // needing its own sync to read back, but only on this stream, so it
    // doesn't block NonBonded or the other bonded/special terms running
    // concurrently on their own streams.
    cudaStream_t m_stream = 0;
  };

} // namespace interaction
