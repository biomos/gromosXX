/**
 * @file read_input.cc
 * implementation of function read_input
 */

#include <util/stdheader.h>
#include <fstream>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

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


int io::read_input(io::Argument const & args,
		   topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   algorithm::Algorithm_Sequence & md_seq)
{
  // create an in_parameter

  // double start = util::now();

  std::ifstream input_file, topo_file, pttopo_file, conf_file;
  
  try{
    input_file.open(args["input"].c_str());
  }
  catch(std::string s){
    io::messages.add("opening input failed", "read_input",
		     io::message::critical);
    return -1;
  }
  if (!input_file.is_open()){
    std::cout << "\n\ncould not open " << args["input"] << "!\n" << std::endl;
    io::messages.add("opening input failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("parameter read from " + args["input"],
		   "read input",
		   io::message::notice);
  
  io::In_Parameter ip(input_file);
  ip.read(sim.param());
  sim.time_step_size() = sim.param().step.dt;
  sim.time() = sim.param().step.t0;

  // std::cout << "timing (parameter) : "
  // << util::now() - start << std::endl;
  // start = util::now();
  
  try{
    topo_file.open(args["topo"].c_str());
  }
  catch(std::string s){
    io::messages.add("opening topology failed", "read_input",
		     io::message::critical);
    return -1;
  }
  if (!topo_file.is_open()){
    std::cout << "\n\ncould not open " << args["topo"] << "!\n" << std::endl;
    io::messages.add("opening topology failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("topology read from " + args["topo"],
		   "read input",
		   io::message::notice);
  
  io::In_Topology it(topo_file);
  it.read(topo, sim.param());

  // std::cout << "timing (topology) : "
  // << util::now() - start << std::endl;
  // start = util::now();

  if(sim.param().perturbation.perturbation){
    if(args.count("pttopo")<1){
      io::messages.add("No perturbation topology specified",
		       "read_input", io::message::critical);
      throw std::string("No perturbation topology specified");
    }
    try{
      pttopo_file.open(args["pttopo"].c_str());
    }
    catch(std::string s){
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::critical);
      return -1;
    }
    if (!pttopo_file.is_open()){
      std::cout << "\n\ncould not open " << args["pttopo"] << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::critical);
      return -1;
    }
    
    io::messages.add("perturbation topology read from " + args["pttopo"],
		     "read input",
		     io::message::notice);
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.read(topo, sim.param());

    // std::cout << "timing (perturbation topology) : "
    // << util::now() - start << std::endl;
    // start = util::now();

  }

  // do this after reading in a perturbation topology
  sim.multibath().calculate_degrees_of_freedom(topo);

  // std::cout << "timing (calculate dof) : "
  // << util::now() - start << std::endl;
  // start = util::now();

  try{
    conf_file.open(args["conf"].c_str());
  }
  catch(std::string s){
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
    return -1;
  }
  if (!conf_file.is_open()){
    std::cout << "\n\ncould not open " << args["conf"] << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("configuration read from " + args["conf"],
		   "read input",
		   io::message::notice);
  
  io::In_Configuration ic(conf_file);
  ic.read(conf, topo, sim.param());

  // std::cout << "timing (read configuration) : "
  // << util::now() - start << std::endl;
  // start = util::now();
  
  // and create the algorithms
  // (among them the forcefield!)
  algorithm::create_md_sequence(md_seq, topo, conf, sim, it);

  // std::cout << "timing (create algorithm) : "
  // << util::now() - start << std::endl;
  
  return 0;
}
