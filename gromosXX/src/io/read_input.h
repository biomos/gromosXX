/**
 * @file read_input.h
 * input parameters
 */

#ifndef INCLUDED_READ_INPUT_H
#define INCLUDED_READ_INPUT_H

namespace io
{
  /**
   * read in data and setup an md simulation.
   */
  int read_input(io::Argument const &args,
		 topology::Topology &topo,
		 configuration::Configuration &conf,
		 simulation::Simulation & sim,
		 algorithm::Algorithm_Sequence &md_seq);
  
}

#endif
