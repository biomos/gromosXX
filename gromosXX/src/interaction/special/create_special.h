/**
 * @file create_special.h
 * create the bonded terms.
 */

#ifndef INCLUDED_CREATE_SPECIAL_H
#define INCLUDED_CREATE_SPECIAL_H

namespace interaction
{
  int create_special(interaction::Forcefield & ff,
		     topology::Topology const & topo,
		     simulation::Parameter const & param,
		     std::ostream & os = std::cout,
		     bool quiet = false);

}

#endif
