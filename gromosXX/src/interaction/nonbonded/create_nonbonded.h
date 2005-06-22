/**
 * @file create_nonbonded.h
 */

#ifndef INCLUDED_CREATE_NONBONDED_H
#define INCLUDED_CREATE_NONBONDED_H

namespace topology
{
  class Topology;
}
namespace simulation
{
  class Simulation;
}
namespace io
{
  class IFP;
}

namespace interaction
{
  class Forcefield;
  
  /**
   * create the nobonded terms for a
   * Gromos96 like simulation.
   */
  int create_g96_nonbonded(interaction::Forcefield & ff,
			   topology::Topology const & topo,
			   simulation::Simulation const & sim,
			   io::IFP & it,
			   std::ostream & os = std::cout,
			   bool quiet = false);

} // interaction

#endif

