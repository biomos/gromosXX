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
 * @file nosehoover_thermostat_cpu.cc
 * methods of the Nose-Hoover Thermostat -- CPU backend.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "nosehoover_thermostat.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

template<>
int algorithm::NoseHoover_Thermostat<util::cpuBackend>
::apply
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  m_timer.start(sim);

  assert(sim.param().multibath.algorithm > 0);

  if (sim.param().multibath.algorithm == 1)
    calc_scaling(topo, conf, sim);
  else
    calc_chain_scaling(topo, conf, sim);

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  scale(topo, conf, sim);

  m_timer.stop();

  return 0;

}

// explicit instantiation for linker
template class algorithm::NoseHoover_Thermostat<util::cpuBackend>;
