/**
 * @file read_special.cc
 * implementation of function read_special
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
#include <io/configuration/inframe.h>
#include <io/configuration/in_configuration.h>
#include <io/topology/in_topology.h>
#include <io/topology/in_perturbation.h>
#include <io/parameter/in_parameter.h>
#include <io/topology/in_posres.h>

#include "read_special.h"

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE read_input


int io::read_special(io::Argument const & args,
		     topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim)
{

  // POSRES
  if (sim.param().posrest.posrest){
    std::ifstream posres_file;
  
    try{
      posres_file.open(args["posres"].c_str());
    }
    catch(std::string s){
      s = "opening posres file failed!\n" + s;
      throw s;
    }
    io::messages.add("position restraints read from " + args["posres"],
		     "read special",
		     io::message::notice);
  
    io::In_Posres ip(posres_file);
    ip.read(topo, conf, sim);
    
  } // POSRES
  
  return 0;
}
