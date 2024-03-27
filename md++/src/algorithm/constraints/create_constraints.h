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
 * @file create_constraints.h
 * create the constraints classes necessary
 * for an md run.
 */

#ifndef INCLUDED_CREATE_CONSTRAINTS_H
#define INCLUDED_CREATE_CONSTRAINTS_H

namespace topology
{
  class Topology;
}
namespace Simulation
{
  class Simulation;
}
namespace io
{
  class In_Topology;
}

namespace algorithm
{
  class Algorithm_Sequence;
  
  int create_constraints(algorithm::Algorithm_Sequence & md_seq,
			 topology::Topology & topo,
			 simulation::Simulation & sim,
			 io::In_Topology &it,
			 bool quiet = false);
  
}

#endif
