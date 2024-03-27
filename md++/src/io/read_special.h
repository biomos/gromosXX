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
 * @file read_special.h
 * read special files
 */

#ifndef INCLUDED_READ_SPECIAL_H
#define INCLUDED_READ_SPECIAL_H

namespace io
{
  int read_special(io::Argument const &args,
		   topology::Topology &topo,
		   configuration::Configuration &conf,
		   simulation::Simulation & sim,
		   std::ostream & os = std::cout,
		   bool quiet = false);
  
}

#endif
