/**
 * @file read_input.h
 * input parameters
 */

#ifndef INCLUDED_READ_INPUT_H
#define INCLUDED_READ_INPUT_H

namespace io
{
  int read_input(io::Argument const &args,
		 simulation::Parameter &param,
		 topology::Topology &topo,
		 configuration::Configuration &conf,
		 algorithm::Algorithm_Sequence &md_seq);
  
}

#endif
