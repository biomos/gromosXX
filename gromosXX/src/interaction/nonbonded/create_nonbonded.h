/**
 * @file create_nonbonded.h
 */

#ifndef INCLUDED_CREATE_NONBONDED_H
#define INCLUDED_CREATE_NONBONDED_H

namespace interaction
{

  void create_g96_nonbonded(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Parameter const & param,
			    io::In_Topology & it);

} // interaction

#endif

