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
 * @file analyze.cc
 * contains the implementation
 * for trajectory analysis
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../io/blockinput.h"
#include "../../io/instream.h"
#include "../../io/configuration/inframe.h"
#include "../../io/configuration/in_configuration.h"

#include "../../math/periodicity.h"
#include "../../util/template_split.h"

#include "analyze.h"

#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

algorithm::Analyze_Step::Analyze_Step(std::string trajectory)
  : Algorithm("AnalyzeStep")
{
  DEBUG(7, "reading configuration");
  m_trajectory_file.open(trajectory.c_str());

  if (!m_trajectory_file.is_open()){
    std::cout << "\n\ncould not open " << trajectory << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
  }

  io::messages.add("analyzation trajectory read from " + trajectory,
		   "read input",
		   io::message::notice);
  
  // and associate
  m_trajectory.stream(m_trajectory_file);

}



template<math::boundary_enum b>
void _put_cg_into_box(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  math::Periodicity<b> p(conf.current().box);
  // just call the function from the periodicity
  p.put_chargegroups_into_box(conf, topo);
}


/**
 * analyzation step init
 */
int algorithm::Analyze_Step
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet) {
    os << "TRAJECTORY ANALYSIS\n"
       << "\tCoordinates will be read in from a previously generated trajectory\n"
       << "\tinstead of being calculated.\n"
       << "\n"
       << "\tNTSTR has to be set to the NTWX parameter\n"
       << "\twhich was used to produce the trajectory being analyzed!\n"
       << "\tmake sure NTWE, NTWX, NTWG, NTWB and NTWF are a multiple\n"
       << "\tof or equal to NTSTR (or 0)!\n";
    if (sim.param().analyze.no_constraints) {
      os << "\n\tConstraints turned off by NTSHK=2.\n";
    }
    os<< "\t\n";
  }
    
 
  // save initial configuration as old() state
  conf.exchange_state();

  if (m_trajectory.read_next(topo, conf, sim, os, false)){

    if (sim.param().analyze.copy_pos){
      os << "\tcopying current position to old\n"
       << "\tvelocities are set to 0.\n";
      DEBUG(9, "analyze: copy current() pos to old()");
      conf.old().pos = conf.current().pos;
    }

    conf.old().vel = 0;
    conf.current().vel = 0;

    SPLIT_BOUNDARY(_put_cg_into_box, topo, conf, sim);    
  
    if (!quiet)  os << "END\n";
    
    return 0;
  }


  return E_UNSPECIFIED;
}

/**
 * trajectory analyzation step.
 */
int algorithm::Analyze_Step
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  // just overwrite current conf
  DEBUG(8, "analyze: reading next frame!");
  conf.exchange_state();

  if (unsigned(sim.param().step.number_of_steps) == sim.steps() + sim.param().analyze.stride){
    // last frame, there is no more frame to read from the trajectory, but we still need to
    // do the step to write out the energies
    return 0;
  } else if (m_trajectory.read_next(topo, conf, sim, std::cout, false)){

    if (sim.param().analyze.copy_pos){
      DEBUG(9, "analyze: copy current() pos to old()");
      conf.old().pos = conf.current().pos;
    }

    conf.old().vel = 0;
    conf.current().vel = 0;
    SPLIT_BOUNDARY(_put_cg_into_box, topo, conf, sim);
    
    return 0;
  } else  {
      io::messages.add("no more frames in the analyzed trajectory",
		   "read input",
		   io::message::error);
    return E_UNSPECIFIED;
  }
}

