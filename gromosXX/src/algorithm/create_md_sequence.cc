/**
 * @file create_md_sequence.cc
 */

#include <util/stdheader.h>
#include <fstream>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <io/argument.h>
#include <io/blockinput.h>
#include <io/instream.h>
#include <io/topology/in_topology.h>

#include <algorithm/algorithm.h>
#include <algorithm/algorithm_sequence.h>
#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen.h>

#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <math/periodicity.h>
#include <algorithm/constraints/shake.h>

#include "create_md_sequence.h"


#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

int algorithm::create_md_sequence(algorithm::Algorithm_Sequence &md_seq,
				  topology::Topology &topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  io::In_Topology &it)
{

  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
  interaction::create_g96_forcefield(*ff, topo, sim.param(), it);

  // construct the md algorithm
  md_seq.push_back(ff);
  md_seq.push_back(new algorithm::Leap_Frog_Velocity);
  md_seq.push_back(new algorithm::Leap_Frog_Position);

  // SHAKE
  DEBUG(7, "SHAKE?");
  if (sim.param().system.nsm || sim.param().shake.ntc > 1){
    DEBUG(8, "\tyes, we need it");
    switch(sim.param().pcouple.virial){
      case math::no_virial:
      case math::molecular_virial:
	{
	  DEBUG(8, "\twith no virial");
	  algorithm::Shake<math::no_virial> * s = 
	    new algorithm::Shake<math::no_virial>
	    (sim.param().shake.tolerance);
	  it.read_harmonic_bonds(s->parameter());
	  md_seq.push_back(s);
	  break;
	}
      case math::atomic_virial:
	DEBUG(8, "\twith atomic virial");
	  algorithm::Shake<math::atomic_virial> * s = 
	    new algorithm::Shake<math::atomic_virial>
	    (sim.param().shake.tolerance);
	  it.read_harmonic_bonds(s->parameter());
	  md_seq.push_back(s);
	break;
    }
  }

  // temperature calculation
  {
    algorithm::Temperature_Calculation * tcalc =
      new algorithm::Temperature_Calculation;
    // calculate initial temperature
    tcalc->apply(topo, conf, sim);

    md_seq.push_back(tcalc);
    
  }
  
  return 0;

}

