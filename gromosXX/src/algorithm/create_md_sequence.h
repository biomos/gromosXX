/**
 * @file create_md_sequence.h
 * create the algorithm and the forcefield
 * for an md run.
 */

#ifndef INCLUDED_CREATE_MD_SEQUENCE_H
#define INCLUDED_CREATE_MD_SEQUENCE_H

namespace algorithm
{
  int create_md_sequence(algorithm::Algorithm_Sequence &md_seq,
			 topology::Topology &topo,
			 simulation::Parameter &param,
			 io::In_Topology &it);
}

#endif

  


