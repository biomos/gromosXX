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
#include <io/topology/in_distanceres.h>
#include <io/topology/in_dihrest.h>
#include <io/topology/in_jvalue.h>
#include <io/topology/in_friction.h>

#include "read_special.h"

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE read_input

/**
 * WARNING: configuration has not been read yet
 */
int io::read_special(io::Argument const & args,
		     topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os,
		     bool quiet)
{
  // POSRES
  if (sim.param().posrest.posrest){
    io::igzstream posresspec_file;
  
    if (args.count("posresspec") != 1){
      io::messages.add("position restraints: no data file specified (use @posresspec)",
		       "read special", io::message::error);
    } else {
      posresspec_file.open(args["posresspec"].c_str());
      if (!posresspec_file.is_open()){
	io::messages.add("opening posresspec file failed!\n",
			 "read_special", 
			 io::message::error);
      } else {
        io::messages.add("position restraints specifciation read from " + args["posresspec"],
                "read special", io::message::notice);

        io::In_Posresspec ip(posresspec_file);
        ip.quiet = quiet;

        ip.read(topo, sim, os);
      }
    }
    
    io::igzstream posres_file;

    // check whether we also need the position restraints file containing the
    // positions and B-factors
    if (sim.param().posrest.posrest == simulation::posrest_bfactor ||
        sim.param().posrest.read) {

      if (args.count("posres") != 1) {
        io::messages.add("position restraints: no data file specified (use @posres)",
                "read special", io::message::error);
      } else {
        posres_file.open(args["posres"].c_str());
        if (!posres_file.is_open()) {
          io::messages.add("opening posresspec file failed!\n",
                  "read_special",
                  io::message::error);
        } else {
          io::messages.add("position restraints specifciation read from " + args["posres"],
                  "read special", io::message::notice);

          io::In_Posres ip(posres_file);
          ip.quiet = quiet;

          ip.read(topo, sim, os);
        }
      }
    }
  } // POSRES

  // DISTANCERES
  if (sim.param().distanceres.distanceres){
    io::igzstream distanceres_file;

    if (args.count("distrest") != 1){
      io::messages.add("distance restraints: no data file specified (use @distrest)",
		       "read special", io::message::error);
    } else {
      distanceres_file.open(args["distrest"].c_str());
      if (!distanceres_file.is_open()){
	io::messages.add("opening distanceres file failed!\n",
			 "read_special", io::message::error);
      } else {
        io::messages.add("distance restraints read from " + args["distrest"],
                "read special", io::message::notice);

        io::In_Distanceres ip(distanceres_file);
        ip.quiet = quiet;

        ip.read(topo, sim, os);
      }
    }    
  } // DISTANCERES

  // DIHREST
  if (sim.param().dihrest.dihrest){
    io::igzstream dihrest_file;

    if (args.count("dihrest") != 1){
      io::messages.add("dihedral restraints: no data file specified (use @dihrest)",
		       "read special", io::message::error);
    } else{
      dihrest_file.open(args["dihrest"].c_str());
      if (!dihrest_file.is_open()){
	io::messages.add("opening dihrest file '" + args["dihrest"] + "'failed!\n",
			 "read_special", io::message::error);
      } else {
        io::messages.add("dihedral restraints read from " + args["dihrest"],
                "read special", io::message::notice);

        io::In_Dihrest ip(dihrest_file);
        ip.quiet = quiet;

        ip.read(topo, sim, os);
      }
    }    
  } // DIHREST

  // J-Value restraints
  if (sim.param().jvalue.mode != simulation::restr_off){
    io::igzstream jval_file;
    
    if (args.count("jval") != 1){
      io::messages.add("jvalue restraints: no data file specified (use @jval)",
		       "read special", io::message::error);
    } else {
      jval_file.open(args["jval"].c_str());
      if (!jval_file.is_open()){
	io::messages.add("opening jvalue restraints file failed!\n",
			 "read_special", io::message::error);
      } else {
        io::messages.add("jvalue restraints read from " + args["jval"],
                "read special", io::message::notice);

        io::In_Jvalue ij(jval_file);
        ij.quiet = quiet;

        ij.read(topo, conf, sim, os);
      }
    }
  } // JVALUE
  
    // FRICTION
  if (sim.param().stochastic.sd && sim.param().stochastic.ntfr == 2){
    io::igzstream friction_file;
  
    if (args.count("friction") != 1){
      io::messages.add("friction specification: no data file specified (use @friction)",
		       "read special", io::message::error);
    }
    else{
      friction_file.open(args["friction"].c_str());
      if (!friction_file.is_open()) {
        io::messages.add("opening friction file failed!\n",
                "read_special", io::message::error);
      } else {
        io::messages.add("atomic friction coefficients read from " + args["friction"],
                "read special", io::message::notice);

        io::In_Friction infr(friction_file);
        infr.quiet = quiet;

        infr.read(topo, sim, os);
      }
    }
  } // FRICTION

  // check for errors and abort
  if (io::messages.contains(io::message::error) ||
      io::messages.contains(io::message::critical))
    return -1;
  
  return 0;
}
