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
 * @file nosehoover_thermostat_gpu.cc
 * Nose-Hoover thermostat -- GPU backend. Plain C++, no kernel syntax;
 * reuses `gpu::apply_thermostat_velocity_scale()`
 * (thermostat_velocity_scale.h) verbatim -- the same O(num_atoms)
 * velocity-scaling helper Berendsen_Thermostat<gpuBackend> uses, since
 * both thermostats inherit the exact same `Thermostat::scale()`
 * formula. Only calc_scaling()/calc_chain_scaling() (O(num_baths)
 * scalar math, shared/backend-agnostic, header-inline) differ between
 * them.
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/algorithm/temperature/thermostat_velocity_scale.h"

#include "nosehoover_thermostat.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE temperature

template<>
int algorithm::NoseHoover_Thermostat<util::gpuBackend>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start(sim);

  assert(sim.param().multibath.algorithm > 0);

  if (sim.param().multibath.algorithm == 1)
    calc_scaling(topo, conf, sim);
  else
    calc_chain_scaling(topo, conf, sim);

  gpu::apply_thermostat_velocity_scale(topo, conf, sim);

  m_timer.stop();
  return 0;
}

// explicit instantiation for linker
template class algorithm::NoseHoover_Thermostat<util::gpuBackend>;
