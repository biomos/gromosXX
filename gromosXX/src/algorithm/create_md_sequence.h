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
namespace configuration{
	class Configuration;
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
	 */
  int create_md_sequence(algorithm::Algorithm_Sequence & md_seq,
			 topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 io::In_Topology &it,
			 std::ostream & os = std::cout);
}

#endif

  


