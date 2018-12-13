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
