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
 * @file temperature_calculation_gpu.cc
 * calculates the temperature -- GPU backend. Plain C++, no kernel
 * syntax (groalgorithm has no CUDA language enabled -- same
 * constraint leap_frog_gpu.cc/remove_com_motion_gpu.cc already work
 * within); the real __global__ kernels live in
 * gpu/cuda/algorithm/temperature/temperature_kernels.cu, compiled into
 * grocuda.
 *
 * Mirrors temperature_calculation_cpu.cc's per-molecule loop exactly,
 * just fed by GPU-reduced per-temperature-group sums (mass, momentum,
 * self-energy, for both conf.current().vel and conf.old().vel) instead
 * of a CPU loop over atoms per molecule -- see temperature_kernels.h's
 * doc comment for the algebra that makes two calls to the same
 * reduction kernel sufficient (no separate "averaged" reduction
 * needed).
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/temperature/temperature_kernels.h"

#include "temperature_calculation.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE temperature

namespace {

  /**
   * Per-atom temperature-group index, built once from
   * topo.temperature_groups(). NOT the same CSR convention as
   * topo.energy_groups()/atom_energy_group(): temperature_groups() has
   * num_groups+1 entries with a leading 0 and EXCLUSIVE end boundaries
   * (confirmed directly from in_topology.cc's TEMPERATUREGROUPS parsing
   * -- "if size()==0 { push_back(0); push_back(num_solute_atoms); }",
   * num_solute_temperature_groups() = size()-1) -- group g spans atoms
   * [tg[g], tg[g+1]). Function-local static: topology is static for a
   * normal run (same assumption as CUDA_Pairlist_Algorithm_Impl's
   * m_iac/m_charge/m_atom_energy_group).
   */
  const gpu::cuvector<unsigned> & atom_temp_group(const topology::Topology & topo) {
    static gpu::cuvector<unsigned> group_index;
    static unsigned cached_num_atoms = 0;
    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    if (cached_num_atoms == num_atoms) return group_index;

    group_index.resize(num_atoms);
    const std::vector<unsigned int> & tg = topo.temperature_groups();
    unsigned group = 0;
    for (unsigned atom = 0; atom < num_atoms; ++atom) {
      if (group + 1 < tg.size() && atom == tg[group + 1]) ++group;
      group_index[atom] = group;
    }
    cached_num_atoms = num_atoms;
    return group_index;
  }

  unsigned num_temperature_groups(const topology::Topology & topo) {
    return static_cast<unsigned>(topo.temperature_groups().size()) - 1;
  }

}

template<>
int algorithm::Temperature_Calculation<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start(sim);

  if (sim.param().perturbation.perturbation) {
    io::messages.add(
        "Temperature_Calculation<gpuBackend>: perturbed kinetic energy "
        "derivatives are not supported (out of scope, matching every "
        "other perturbation gate in the CUDA path).",
        "Temperature_Calculation", io::message::error);
    m_timer.stop();
    return 1;
  }

  // zero previous (temperature scaling) energies
  conf.old().energies.zero(false, true);
  for (unsigned i = 0; i < unsigned(sim.multibath().size()); ++i)
    sim.multibath().bath(i).ekin = 0.0;

  const unsigned num_atoms  = static_cast<unsigned>(topo.num_atoms());
  const unsigned num_groups = num_temperature_groups(topo);
  const gpu::cuvector<unsigned> & group_index = atom_temp_group(topo);

  // Own stream: both reductions (new_sums from current().vel, old_sums
  // from old().vel) queue on it back-to-back with no CPU involvement
  // between them, then one stream-scoped (not device-wide) wait below.
  static cudaStream_t stream = 0;
  if (stream == 0) cudaStreamCreate(&stream);

  // MIRROR_VEL freshness covers current()+old() together (mark_gpu_dirty
  // always moves both, matching copy_pos_vel_to_device()'s own
  // granularity) -- now that CudaManager::exchange_mirror_state() keeps
  // the GPU mirror's current/old halves in lockstep with conf.
  // exchange_state() (fixed alongside this), view.old().vel here is
  // already correct GPU-resident data: no algorithm swaps state again
  // between Leap_Frog_Velocity and this call, so it's exactly what CPU's
  // conf.old().vel holds too. This replaces what used to be a per-atom
  // host upload loop (old_vel_gpu) every single call.
  gpu::Configuration::View conf_view =
      sim.cuda().configuration_view(conf, gpu::MIRROR_VEL, stream);
  const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

  static gpu::cuvector<double> new_sums, old_sums;
  if (new_sums.size() < 5u * num_groups) new_sums.resize(5u * num_groups);
  if (old_sums.size() < 5u * num_groups) old_sums.resize(5u * num_groups);

  gpu::launch_group_velocity_reduce(conf_view.current().vel, topo_view.mass,
                                     group_index.data(), num_atoms, num_groups,
                                     new_sums.data(), stream);
  gpu::launch_group_velocity_reduce(conf_view.old().vel, topo_view.mass,
                                     group_index.data(), num_atoms, num_groups,
                                     old_sums.data(), stream);
  cudaStreamSynchronize(stream);

  unsigned ir_bath = 0, com_bath = 0;

  for (unsigned g = 0; g < num_groups; ++g) {
    const double mass = new_sums[5u*g + 0];
    if (mass <= 0.0) continue; // empty group (shouldn't happen, but stay safe)

    const math::Vec new_com_v(new_sums[5u*g+1]/mass, new_sums[5u*g+2]/mass, new_sums[5u*g+3]/mass);
    const math::Vec old_com_v(old_sums[5u*g+1]/mass, old_sums[5u*g+2]/mass, old_sums[5u*g+3]/mass);
    const math::Vec com_v_avg = 0.5 * (new_com_v + old_com_v);

    const double new_com_ekin = 0.5 * mass * abs2(new_com_v);
    const double new_ekin     = 0.5 * new_sums[5u*g+4];
    const double com_ekin_avg = 0.5 * mass * abs2(com_v_avg);
    const double ekin_avg     = 0.25 * (new_sums[5u*g+4] + old_sums[5u*g+4]);

    // Group g spans atoms [tg[g], tg[g+1]) -- tg[g] is its first atom,
    // which identifies the group for in_bath() (same convention as the
    // CPU loop's *tg_it.begin()).
    const std::vector<unsigned int> & tg = topo.temperature_groups();
    sim.multibath().in_bath(tg[g], com_bath, ir_bath);

    sim.multibath().bath(com_bath).ekin += new_com_ekin;
    sim.multibath().bath(ir_bath).ekin  += new_ekin - new_com_ekin;

    conf.old().energies.com_kinetic_energy[com_bath] += com_ekin_avg;
    conf.old().energies.ir_kinetic_energy[ir_bath]    += ekin_avg - com_ekin_avg;
  }

  for (size_t i = 0; i < conf.old().energies.kinetic_energy.size(); ++i)
    conf.old().energies.kinetic_energy[i] =
        conf.old().energies.com_kinetic_energy[i] +
        conf.old().energies.ir_kinetic_energy[i];

  m_timer.stop();
  return 0;
}

// explicit instantiation for linker
template class algorithm::Temperature_Calculation<util::gpuBackend>;
