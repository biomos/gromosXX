/**
 * @file read_input.cc
 * implementation of function read_input
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

#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/create_md_sequence.h>

#include <interaction/forcefield/forcefield.h>

#include "read_input.h"
#include "read_special.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE parameter


int io::read_input(io::Argument const & args,
		   topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   algorithm::Algorithm_Sequence & md_seq,
		   std::ostream & os)
{
  // create an in_parameter

  // double start = util::now();

  std::ifstream input_file, topo_file, pttopo_file, conf_file;
  
  input_file.open(args["input"].c_str());

  if (!input_file.is_open()){
    os << "\n\ncould not open " << args["input"] << "!\n" << std::endl;
    io::messages.add("opening input failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("parameter read from " + args["input"],
		   "read input",
		   io::message::notice);
  
  io::In_Parameter ip(input_file);
  ip.read(sim.param(), os);
  sim.time_step_size() = sim.param().step.dt;
  sim.time() = sim.param().step.t0;

  if (sim.param().analyze.analyze){
    if (args.count("anatrj") < 1){
      os << "\n\nno analyzation trajectory specified (@anatrj)\n";
      io::messages.add("\n\nno analyzation trajectory specified (@anatrj)\n",
		       "read_input",
		       io::message::critical);
      sim.param().analyze.analyze = false;
    }
    else
      sim.param().analyze.trajectory = args["anatrj"];
  }

  topo_file.open(args["topo"].c_str());
  if (!topo_file.is_open()){
    os << "\n\ncould not open " << args["topo"] << "!\n" << std::endl;
    io::messages.add("opening topology failed", "read_input",
		     io::message::critical);
  }

  io::messages.add("topology read from " + args["topo"],
		   "read input",
		   io::message::notice);
  
  io::In_Topology it(topo_file);
  it.read(topo, sim.param(), os);

  if(sim.param().perturbation.perturbation){
    if(args.count("pttopo")<1){
      io::messages.add("No perturbation topology specified",
		       "read_input", io::message::critical);
      return -1;
    }

    pttopo_file.open(args["pttopo"].c_str());

    if (!pttopo_file.is_open()){
      os << "\n\ncould not open " << args["pttopo"] << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::critical);
    }
    
    io::messages.add("perturbation topology read from " + args["pttopo"],
		     "read input",
		     io::message::notice);
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.read(topo, sim.param());

  }

  topo.init(sim, os);

  // do this after reading in a perturbation topology
  sim.multibath().calculate_degrees_of_freedom(topo, sim.param().rottrans.rottrans);

  // read in the special data
  DEBUG(7, "reading special data");
  io::read_special(args, topo, conf, sim);

  DEBUG(7, "reading configuration");
  conf_file.open(args["conf"].c_str());

  if (!conf_file.is_open()){
    os << "\n\ncould not open " << args["conf"] << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
    return -1;
  }


  io::messages.add("configuration read from " + args["conf"],
		   "read input",
		   io::message::notice);
  
  io::In_Configuration ic(conf_file);
  ic.read(conf, topo, sim, os);
  // and initialise
  conf.initialise(topo, sim.param());
    
  // and create the algorithms
  // (among them the forcefield!)
  algorithm::create_md_sequence(md_seq, topo, conf, sim, it, os);

  return 0;
}
