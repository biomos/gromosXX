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
 * @file lattice_shift.cc
 * implementation of the lattice shift tracking
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../math/periodicity.h"
#include "../../util/template_split.h"

#include "lattice_shift.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

int algorithm::Lattice_Shift_Tracker::
init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {
  
  // initialize the lattice shifts if not read from configuration
  if (!sim.param().start.read_lattice_shifts) {
    conf.special().lattice_shifts = 0.0;
  }
  
  if (!quiet) { // write some stuff to the output file
    os << "LATTICESHIFTS" << std::endl
       << "    keeping track of lattice shifts." << std::endl;
    
    if (sim.param().start.read_lattice_shifts)
      os << "    reading initial shifts from configuration.";
    else
      os << "    setting initial shifts to zero.";
    
    os << std::endl
       << "END" << std::endl;
  }
  
  return 0;
}

int algorithm::Lattice_Shift_Tracker::
apply(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  DEBUG(6, "keeping track of lattice shifts");
  m_timer.start(sim);
  SPLIT_BOUNDARY(_apply, topo, conf, sim);
  m_timer.stop();
  return 0;
}

template<math::boundary_enum b>
void algorithm::Lattice_Shift_Tracker::
_apply(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  math::Periodicity<b> p(conf.current().box);
  // just call the function from the periodicity
  p.put_chargegroups_into_box_saving_shifts(conf, topo);
}



