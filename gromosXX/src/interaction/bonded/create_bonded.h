/**
 * @file create_bonded.h
 * create the bonded terms.
 */

#ifndef INCLUDED_CREATE_BONDED_H
#define INCLUDED_CREATE_BONDED_H

namespace interaction
{
  void create_g96_bonded(interaction::Forcefield & ff,
			 topology::Topology const & topo,
			 simulation::Parameter const & param,
			 io::In_Topology & it);

}

#endif
