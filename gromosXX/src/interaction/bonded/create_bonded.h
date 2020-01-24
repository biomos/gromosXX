/**
 * @file create_bonded.h
 * create the bonded terms.
 */

#ifndef INCLUDED_CREATE_BONDED_H
#define INCLUDED_CREATE_BONDED_H

namespace topology
{
  class Topology;
}
namespace simulation
{
  class Parameter;
}
namespace io
{
  class IFP;
}

namespace interaction
{
  class Forcefield;

	/**
	 * create the bonded interaction terms for a Gromos96
	 * like simulation.
	 */
  int create_g96_bonded(interaction::Forcefield & ff,
			topology::Topology const & topo,
			simulation::Parameter const & param,
			io::IFP & it,
			std::ostream & os = std::cout,
			bool quiet = false);

}

#endif
