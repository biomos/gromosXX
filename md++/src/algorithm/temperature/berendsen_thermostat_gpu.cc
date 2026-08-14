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
 * @file berendsen_thermostat_gpu.cc
 * Berendsen thermostat -- GPU backend. Plain C++, no kernel syntax; the
 * real __global__ kernels live in
 * gpu/cuda/algorithm/temperature/temperature_kernels.cu.
 *
 * calc_scaling() (inherited, header-inline, no backend-specific logic)
 * stays exactly as-is -- O(num_baths) scalar math. Only the O(num_atoms)
 * velocity scaling (Thermostat::scale(), CPU-only) is replaced here with
 * the uniform per-atom GPU formula (temperature_kernels.h's doc comment
 * has the derivation).
 *
 * Leaves the scaled velocity resident on the GPU mirror
 * (mark_gpu_dirty(), no sync-back) -- this sits between Leap_Frog_
 * Velocity<gpuBackend> and Leap_Frog_Position<gpuBackend> in
 * create_md_sequence.cc, so the whole Velocity -> Thermostat -> Position
 * chain now runs with zero host round trips when temperature coupling
 * is on, syncing back to the CPU only once, at the very end of Position
 * (its existing sync_configuration_from_device() call).
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/temperature/temperature_kernels.h"

#include "berendsen_thermostat.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE temperature

namespace {

  /**
   * Per-atom temperature-group index and (com_bath, ir_bath) indices,
   * built once from topo.temperature_groups() and
   * sim.multibath().in_bath() (the existing CPU-side lookup, reused
   * rather than re-deriving its range logic). Function-local statics,
   * same rationale as Remove_COM_Motion's/Temperature_Calculation's
   * scratch buffers: topology and multibath structure are static for a
   * normal run.
   */
  struct AtomBathArrays {
    gpu::cuvector<unsigned> group_index;
    gpu::cuvector<unsigned> com_bath;
    gpu::cuvector<unsigned> ir_bath;
    unsigned num_groups = 0;
  };

  const AtomBathArrays & atom_bath_arrays(const topology::Topology & topo,
                                           const simulation::Simulation & sim) {
    static AtomBathArrays arrays;
    static unsigned cached_num_atoms = 0;
    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    if (cached_num_atoms == num_atoms) return arrays;

    arrays.group_index.resize(num_atoms);
    arrays.com_bath.resize(num_atoms);
    arrays.ir_bath.resize(num_atoms);

    // temperature_groups() has num_groups+1 entries with a leading 0
    // and EXCLUSIVE end boundaries (see temperature_calculation_gpu.cc's
    // atom_temp_group() for the confirmation from in_topology.cc) --
    // NOT the same convention as energy_groups(). Group g spans atoms
    // [tg[g], tg[g+1]).
    const std::vector<unsigned int> & tg = topo.temperature_groups();
    unsigned group = 0;
    for (unsigned atom = 0; atom < num_atoms; ++atom) {
      if (group + 1 < tg.size() && atom == tg[group + 1]) ++group;
      arrays.group_index[atom] = group;
      unsigned com = 0, ir = 0;
      sim.multibath().in_bath(atom, com, ir);
      arrays.com_bath[atom] = com;
      arrays.ir_bath[atom]  = ir;
    }
    arrays.num_groups = static_cast<unsigned>(tg.size()) - 1;
    cached_num_atoms = num_atoms;
    return arrays;
  }

}

template<>
int algorithm::Berendsen_Thermostat<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start(sim);

  calc_scaling(topo, conf, sim);

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const AtomBathArrays & arrays = atom_bath_arrays(topo, sim);
  const unsigned num_groups = arrays.num_groups;
  const unsigned num_baths  = static_cast<unsigned>(sim.multibath().size());

  gpu::Configuration::View conf_view =
      sim.cuda().configuration_view(conf, gpu::MIRROR_VEL);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  static gpu::cuvector<double> sums;
  if (sums.size() < 5u * num_groups) sums.resize(5u * num_groups);

  gpu::launch_group_velocity_reduce(conf_view.current().vel, topo_view.mass,
                                     arrays.group_index.data(), num_atoms,
                                     num_groups, sums.data());
  cudaDeviceSynchronize();

  static gpu::cuvector<FPL3_TYPE> com_v_per_group;
  if (com_v_per_group.size() < num_groups) com_v_per_group.resize(num_groups);
  for (unsigned g = 0; g < num_groups; ++g) {
    const double mass = sums[5u*g + 0];
    if (mass > 0.0) {
      com_v_per_group[g] = FPL3_TYPE{
          static_cast<FPL_TYPE>(sums[5u*g+1] / mass),
          static_cast<FPL_TYPE>(sums[5u*g+2] / mass),
          static_cast<FPL_TYPE>(sums[5u*g+3] / mass)};
    } else {
      com_v_per_group[g] = FPL3_TYPE{0, 0, 0};
    }
  }

  static gpu::cuvector<double> bath_scale;
  if (bath_scale.size() < num_baths) bath_scale.resize(num_baths);
  for (unsigned b = 0; b < num_baths; ++b)
    bath_scale[b] = sim.multibath()[b].scale;

  gpu::launch_thermostat_scale_apply(
      conf_view.current().vel, arrays.group_index.data(),
      arrays.com_bath.data(), arrays.ir_bath.data(),
      com_v_per_group.data(), bath_scale.data(), num_atoms);

  // Stay GPU-resident -- no sync back here, see file doc comment.
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VEL);

  m_timer.stop();
  return 0;
}

// explicit instantiation for linker
template class algorithm::Berendsen_Thermostat<util::gpuBackend>;
