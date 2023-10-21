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
 * @file prepare_virial.h
 * prepare virial calculation.
 */

#ifndef INCLUDED_PREPARE_VIRIAL_H
#define INCLUDED_PREPARE_VIRIAL_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace util
{
  /**
   * prepare for the virial calculation.
   */
  void prepare_virial(topology::Topology const & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation const & sim);

  /**
   * recover molecular virial from atomic virial
   */
  void atomic_to_molecular_virial(topology::Topology const & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation const & sim);
  
  /**
   * calculate centre of mass, centre of mass kinetic energy
   * of the (sub) molecules
   */
  void centre_of_mass(topology::Topology const & topo,
		      configuration::Configuration & conf,
		      std::vector<math::Vec> & com_pos,
		      std::vector<math::Matrix> & com_ekin,
          	      simulation::Simulation const & sim);

} // util

#endif
