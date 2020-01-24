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
