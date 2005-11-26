/**
 * @file read_special.cc
 * implementation of function read_special
 */

#include <stdheader.h>
#include <fstream>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
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
#include <io/topology/in_distrest.h>
#include <io/topology/in_dihrest.h>
#include <io/topology/in_jvalue.h>

#include "read_special.h"

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE read_input


int io::read_special(io::Argument const & args,
		     topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os,
		     bool quiet)
{

  // POSRES
  if (sim.param().posrest.posrest){
    std::ifstream posres_file;
  
    if (args.count("posres") != 1){
      io::messages.add("position restraints: no data file specified (use @posres)",
		       "read special",
		       io::message::error);
    }
    else{
      posres_file.open(args["posres"].c_str());
      if (!posres_file){
	io::messages.add("opening posres file failed!\n",
			 "read_special", 
			 io::message::error);
      }
      io::messages.add("position restraints read from " + args["posres"],
		       "read special",
		       io::message::notice);
      
      io::In_Posres ip(posres_file);
      ip.quiet = quiet;
      
      ip.read(topo, conf, sim, os);
    }
  } // POSRES

  // DISTREST
  if (sim.param().distrest.distrest){
    std::ifstream distrest_file;

    if (args.count("distrest") != 1){
      io::messages.add("distance restraints: no data file specified (use @distrest)",
		       "read special",
		       io::message::error);
    }
    else{
  
      distrest_file.open(args["distrest"].c_str());
      if (!distrest_file){
	io::messages.add("opening distrest file failed!\n",
			 "read_special", 
			 io::message::error);
      }
      io::messages.add("distance restraints read from " + args["distrest"],
		       "read special",
		       io::message::notice);
      
      io::In_Distrest ip(distrest_file);
      ip.quiet = quiet;
      
      ip.read(topo, conf, sim, os);
    }    
  } // DISTREST

  // DIHREST
  if (sim.param().dihrest.dihrest){
    std::ifstream dihrest_file;

    if (args.count("dihrest") != 1){
      io::messages.add("dihedral restraints: no data file specified (use @dihrest)",
		       "read special",
		       io::message::error);
    }
    else{
      dihrest_file.open(args["dihrest"].c_str());
      if (!dihrest_file){
	io::messages.add("opening dihrest file '" + args["dihrest"] + "'failed!\n",
			 "read_special", 
			 io::message::error);
      }
      io::messages.add("dihedral restraints read from " + args["dihrest"],
		       "read special",
		       io::message::notice);
      
      io::In_Dihrest ip(dihrest_file);
      ip.quiet = quiet;
      
      ip.read(topo, conf, sim, os);
    }    
  } // DIHREST

  // J-Value restraints
  if (sim.param().jvalue.mode != simulation::restr_off){
    std::ifstream jval_file;
    
    if (args.count("jval") != 1){
      io::messages.add("jvalue restraints: no data file specified (use @jval)",
		       "read special",
		       io::message::error);
    }
    else{
      jval_file.open(args["jval"].c_str());
      if (!jval_file){
	io::messages.add("opening jvalue restraints file failed!\n",
			 "read_special", io::message::error);
	return 1;
      }
      io::messages.add("jvalue restraints read from " + args["jval"],
		       "read special",
		       io::message::notice);
      
      io::In_Jvalue ij(jval_file);
      ij.quiet = quiet;
      
      ij.read(topo, conf, sim, os);
    }
  } // JVALUE
  
  return 0;
}
