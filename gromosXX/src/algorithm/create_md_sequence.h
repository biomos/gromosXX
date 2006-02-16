/**
 * @file create_md_sequence.h
 * create the algorithm and the forcefield
 * for an md run.
 */

#ifndef INCLUDED_CREATE_MD_SEQUENCE_H
#define INCLUDED_CREATE_MD_SEQUENCE_H

namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}
namespace io{
	class In_Topology;
}

namespace algorithm
{
	class Algorithm_Sequence;

	/**
	 * create a MD sequence.
	 * creates the necessary algorithm in the correct order and
	 * inserts them into the Algorithm_Sequence md_seq.
	 * @TODO this should be reworked into a
	 * Factory (factory pattern) using sub-factories for the
	 * big thingies like forcefield, maybe algorithms
	 */
  int create_md_sequence(algorithm::Algorithm_Sequence & md_seq,
			 topology::Topology & topo,
			 simulation::Simulation & sim,
			 io::In_Topology &it,
			 std::ostream & os = std::cout,
			 bool quiet = false);
}

#endif

  


