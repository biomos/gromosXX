/**
 * @file create_special.h
 * create the bonded terms.
 */

#ifndef INCLUDED_CREATE_SPECIAL_H
#define INCLUDED_CREATE_SPECIAL_H

namespace interaction
{
  void create_special(interaction::Forcefield & ff,
		      topology::Topology const & topo,
		      simulation::Parameter const & param);

}

#endif
