/**
 * @file read_input.cc
 * implementation of function read_input
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

#include <algorithm/algorithm.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/create_md_sequence.h>

#include <interaction/forcefield/forcefield.h>

#include "read_input.h"

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE read_input


int io::read_input(io::Argument const &args,
		   simulation::Parameter &param,
		   topology::Topology &topo,
		   configuration::Configuration &conf,
		   algorithm::Algorithm_Sequence &md_seq)
{
  // create an in_parameter
  std::ifstream input_file, topo_file, pttopo_file, conf_file;
  
  try{
    input_file.open(args["input"].c_str());
  }
  catch(std::string s){
    s = "opening input failed!\n" + s;
    throw s;
  }
  io::messages.add("parameter read from " + args["input"],
		   "read input",
		   io::message::notice);
  
  io::In_Parameter ip(input_file);
  ip.read(param);
  
  try{
    topo_file.open(args["topo"].c_str());
  }
  catch(std::string s){
    s = "opening topology failed!\n" + s;
    throw s;
  }

  io::messages.add("topology read from " + args["topo"],
		   "read input",
		   io::message::notice);
  
  io::In_Topology it(topo_file);
  it.read(topo, param);

  if(param.perturbation.perturbation){
    if(args.count("pttopo")<1){
      io::messages.add("No perturbation topology specified",
		       "read_input", io::message::critical);
      throw std::string("No perturbation topology specified");
    }
    try{
      pttopo_file.open(args["pttopo"].c_str());
    }
    catch(std::string s){
      s = "opening perturbation topology failed!\n" + s;
    }

    io::messages.add("perturbation topology read from " + args["pttopo"],
		     "read input",
		     io::message::notice);
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.read(topo, param);

  }

  try{
    conf_file.open(args["conf"].c_str());
  }
  catch(std::string s){
    s = "opening configuration failed!\n" + s;
    throw s;
  }

  io::messages.add("configuration read from " + args["conf"],
		   "read input",
		   io::message::notice);
  
  io::In_Configuration ic(conf_file);
  ic.read(conf, topo, param);
  
  // do this after reading in a perturbation topology
  param.multibath.multibath.calculate_degrees_of_freedom(topo);

  // and create the algorithms
  // (among them the forcefield!)
  algorithm::create_md_sequence(md_seq, topo, param, it);
  
  return 0;
}
