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
 * @file cuda_molecular_virial_interaction.h
 * GPU-native molecular-virial correction. Only ever included under
 * USE_CUDA -- see create_forcefield.cc, the only call site (same
 * convention as cuda_quartic_bond_interaction.h/cuda_angle_interaction.h).
 *
 * The CPU reference (Molecular_Virial_Interaction, util/prepare_
 * virial.cc's atomic_to_molecular_virial()) reads conf.current().force
 * and conf.current().virial_tensor as plain CPU arrays and corrects the
 * latter in place. On GPU that CPU-side virial_tensor was never
 * populated with the atomic virial the GPU-native bonded/nonbonded
 * terms accumulate straight into the mirror (a real bug this class
 * replaces the fix for, see git history around "Molecular_Virial_
 * Interaction reads stale/clobbered virial") -- this class instead
 * reads force and applies the correction entirely on-device, into the
 * mirror's own virial_tensor, so nothing ever needs a CPU round trip
 * for this correction and nothing downstream (Pressure_Calculation<
 * gpuBackend>, a later flush) can clobber it.
 *
 * Pressure-group topology (topo.pressure_groups(), a CSR offsets array)
 * is static for a normal run, so the per-atom group_id lookup and group
 * offsets are uploaded once in init(), same convention as every other
 * bonded term's per-atom static lookup tables.
 */

#pragma once

#include "../interaction.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {

  class CUDA_Molecular_Virial_Interaction : public Interaction {
  public:
    CUDA_Molecular_Virial_Interaction() : Interaction("MolecularVirial") {}
    virtual ~CUDA_Molecular_Virial_Interaction();

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int calculate_interactions(topology::Topology & topo,
                                        configuration::Configuration & conf,
                                        simulation::Simulation & sim);

    // Reads force via configuration_view() (GPU-resident, producer-
    // event-ordered), not the plain CPU array -- must NOT trigger
    // Forcefield::calculate_interactions()'s default pre-call FORCE
    // flush (that would force a premature, wasteful publish of whatever
    // an earlier GPU-native term already wrote this step).
    virtual bool needs_fresh_cpu_force() const override { return false; }
    // Writes virial_tensor directly into the mirror via atomicAdd, never
    // touches the CPU array at all -- no CPU-side read/write-back needed
    // (contrast the CPU-only Molecular_Virial_Interaction, which does
    // need both, see needs_fresh_cpu_virial()'s doc comment).
    virtual bool needs_fresh_cpu_virial() const override { return false; }
    virtual bool is_gpu_native() const override { return true; }

  private:
    // CSR pressure-group offsets, size num_groups+1 (topo.pressure_
    // groups() itself, uploaded verbatim -- see molecular_virial_
    // kernels.h's doc comment for why group membership is contiguous).
    gpu::cuvector<unsigned> m_group_offsets;
    // Per-atom group index, size num_atoms -- avoids a per-kernel-call
    // binary search over group_offsets in the second pass.
    gpu::cuvector<unsigned> m_group_id;
    // Per-group mass-weighted unwrapped COM position, size num_groups.
    gpu::cuvector<FPL3_TYPE> m_com_pos;
    // Private virial accumulator (already-negated corrP), row-major
    // b*3+a, matching virial_accumulate_kernels.h's convention -- reused
    // every call, zeroed at the start of each calculate_interactions().
    gpu::cuvector<double> m_virial;

    unsigned m_num_groups = 0;
    unsigned m_num_atoms = 0;
    bool m_initialized = false;

    // Own stream: both kernels launch on this stream, in-order by CUDA's
    // own single-stream ordering guarantee (no explicit sync needed
    // between the group-COM pass and the per-atom correction pass that
    // reads its output) -- and both writes (into m_com_pos, then into
    // the mirror's virial_tensor via virial_accumulate_kernels.h) stay
    // fully GPU-resident, no host round trip anywhere in this class.
    cudaStream_t m_stream = 0;
  };

} // namespace interaction
