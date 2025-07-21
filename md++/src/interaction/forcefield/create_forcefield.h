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
 * @file create_forcefield.h
 */

#ifndef INCLUDED_CREATE_FORCEFIELD_H
#define INCLUDED_CREATE_FORCEFIELD_H

namespace topology
{
  class Topology;
}
namespace simulation
{
  class Simulation;
}
namespace io
{
  class IFP;
}

namespace interaction
{
  /**
   * create a Gromos96 (like) forcefield.
   */
  int create_g96_forcefield(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Simulation & sim,
			    io::IFP & it,
			    std::ostream & os = std::cout,
			    bool quiet = false);
    
  
}

#endif
