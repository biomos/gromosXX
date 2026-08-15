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
 * stays exactly as-is -- O(num_baths) scalar math. The O(num_atoms)
 * velocity scaling is `gpu::apply_thermostat_velocity_scale()`
 * (thermostat_velocity_scale.h), shared with
 * NoseHoover_Thermostat<gpuBackend> since both inherit the exact same
 * `Thermostat::scale()` formula -- see that header's doc comment for
 * the zero-round-trip GPU-mirror-residency reasoning.
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/algorithm/temperature/thermostat_velocity_scale.h"

#include "berendsen_thermostat.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE temperature

template<>
int algorithm::Berendsen_Thermostat<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start(sim);

  calc_scaling(topo, conf, sim);
  gpu::apply_thermostat_velocity_scale(topo, conf, sim);

  m_timer.stop();
  return 0;
}

// explicit instantiation for linker
template class algorithm::Berendsen_Thermostat<util::gpuBackend>;
