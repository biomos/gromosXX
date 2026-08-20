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
 * @file thermostat_velocity_scale.h
 * The O(num_atoms) half of `Thermostat::scale()`, ported to GPU once
 * and shared by every `Thermostat`-derived `gpuBackend` specialization
 * (berendsen_thermostat_gpu.cc, nosehoover_thermostat_gpu.cc): reduce
 * each temperature group's centre-of-mass velocity, then apply the
 * per-atom scale formula. Callers are responsible for the O(num_baths)
 * scalar part (each thermostat's own `calc_scaling()`/
 * `calc_chain_scaling()`, already backend-agnostic and header-inline)
 * *before* calling this.
 *
 * Deliberately leaves the scaled velocity resident on the GPU mirror
 * (`mark_gpu_dirty()`, no sync-back) -- every caller sits between
 * `Leap_Frog_Velocity<gpuBackend>` and `Leap_Frog_Position<gpuBackend>`
 * in `create_md_sequence.cc`, so this keeps that three-algorithm chain
 * running with zero host round trips when temperature coupling is on.
 */

#pragma once

#include "gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/algorithm/temperature/temperature_kernels.h"
#include "gpu/cuda/algorithm/temperature/atom_bath_arrays.h"

namespace gpu {

  inline void apply_thermostat_velocity_scale(topology::Topology & topo,
                                               configuration::Configuration & conf,
                                               simulation::Simulation & sim) {
    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const gpu::AtomBathArrays & arrays = gpu::atom_bath_arrays(topo, sim);
    const unsigned num_groups = arrays.num_groups;
    const unsigned num_baths  = static_cast<unsigned>(sim.multibath().size());

    // Own stream, function-local static (shared by every Thermostat-
    // derived gpuBackend caller through this one function) -- the
    // reduction's result (com_v_per_group, computed below) is a
    // host-side input to launch_thermostat_scale_apply(), an inherent
    // sequential dependency the event mechanism can't route around, but
    // scoping the wait to this stream (not cudaDeviceSynchronize())
    // keeps it from blocking whatever else is running concurrently.
    static cudaStream_t stream = 0;
    if (stream == 0) cudaStreamCreate(&stream);

    gpu::Configuration::View conf_view =
        sim.cuda().configuration_view(conf, gpu::MIRROR_VEL, stream);
    const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

    static gpu::cuvector<double> sums;
    if (sums.size() < 5u * num_groups) sums.resize(5u * num_groups);

    gpu::launch_group_velocity_reduce(conf_view.current().vel, topo_view.mass,
                                       arrays.group_index.data(), num_atoms,
                                       num_groups, sums.data(), stream);
    cudaStreamSynchronize(stream);

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
        com_v_per_group.data(), bath_scale.data(), num_atoms, stream);

    // Stay GPU-resident -- no sync back here, see this file's doc comment.
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VEL, stream);
  }

} // namespace gpu
