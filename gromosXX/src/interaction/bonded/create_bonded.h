/**
 * @file create_bonded.h
 * create the bonded terms.
 */

#ifndef INCLUDED_CREATE_BONDED_H
#define INCLUDED_CREATE_BONDED_H

namespace interaction
{
  int create_g96_bonded(interaction::Forcefield & ff,
			topology::Topology const & topo,
			simulation::Parameter const & param,
			io::IFP & it,
			bool quiet = false);

}

#endif
