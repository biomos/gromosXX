/**
 * @file create_nonbonded.h
 */

#ifndef INCLUDED_CREATE_NONBONDED_H
#define INCLUDED_CREATE_NONBONDED_H

namespace interaction
{

  int create_g96_nonbonded(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   configuration::Configuration const & conf,
			   io::In_Topology & it,
			   bool quiet = false);

} // interaction

#endif

