/**
 * @file create_constraints.h
 * create the constraints classes necessary
 * for an md run.
 */

#ifndef INCLUDED_CREATE_CONSTRAINTS_H
#define INCLUDED_CREATE_CONSTRAINTS_H

namespace algorithm
{
  int create_constraints(algorithm::Algorithm_Sequence & md_seq,
			 topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 io::In_Topology &it);
}

#endif
