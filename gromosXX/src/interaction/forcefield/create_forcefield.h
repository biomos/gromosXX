/**
 * @file create_forcefield.h
 */

#ifndef INCLUDED_CREATE_FORCEFIELD_H
#define INCLUDED_CREATE_FORCEFIELD_H

namespace io
{
  class IFP;
}

namespace interaction
{
  
  /**
   * create a Gromos96 (like) forcefield.
   */
  int create_g96_forcefield(interaction::Forcefield & ff,
			    topology::Topology const & topo,
			    simulation::Simulation const & sim,
			    configuration::Configuration const & conf,
			    io::IFP & it,
			    bool quiet = false);
    
  
}

#endif
