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
 * @file create_simulation.h
 * create a minimum simulation
 */

#ifndef INCLUDED_CREATE_SIMULATION_H
#define INCLUDED_CREATE_SIMULATION_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}
namespace algorithm{
	class Algorithm_Sequence;
}
namespace io{
	class In_Topology;
}

namespace util
{

  /**
   * contains all important entities
   * for an md simulation.
   */
  struct simulation_struct
  {
    topology::Topology topo;
    configuration::Configuration conf;
    simulation::Simulation sim;
    algorithm::Algorithm_Sequence md;
  };
  
  /**
   * create a minimum simulation environment.
   * intended for tests (and maybe Gromos++)
   */
  int create_simulation(std::string topo,
			std::string pttopo,
			std::string conf,
			std::string param,
			util::simulation_struct & sim,
			io::In_Topology & in_topo,
			std::string distanceres = "",
			std::string angrest = "",
			std::string dihrest = "",
                        std::string xray = "",
                        std::string qmmm = "",
			std::string lud = "",
			std::string led = "",
                        std::string order = "",
			bool quiet = false);
}

#endif
