/**
 * @file create_simulation.cc
 * implementation of function create_simulation
 * to easily create a simulation from (non-complete) data.
 */

#include <util/stdheader.h>
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

#include "create_simulation.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

int util::create_simulation(std::string topo,
			    std::string pttopo,
			    std::string conf,
			    std::string param,
			    util::simulation_struct & sim,
			    io::In_Topology & in_topo,
			    bool quiet)
{

  // a topology is needed
  if (topo == ""){
    io::messages.add("No topology specified",
		     "create_simulation",
		     io::message::error);
    return 1;
  }

  std::ifstream input_file, topo_file, pttopo_file, conf_file;
  
  // if we got a parameter file, try to read it...
  if (param != ""){
    try{
      input_file.open(param.c_str());
    }
    catch(std::string s){
      io::messages.add("opening input failed", "read_input",
		       io::message::error);
      return -1;
    }
    if (!input_file.is_open()){
      std::cout << "\n\ncould not open " << param << "!\n" << std::endl;
      io::messages.add("opening input failed", "read_input",
		       io::message::error);
      return -1;
    }
    
    io::In_Parameter ip(input_file);
    ip.quiet = quiet;
    ip.read(sim.sim.param());

    sim.sim.time_step_size() = sim.sim.param().step.dt;
    sim.sim.time() = sim.sim.param().step.t0;

  }
  
  
  try{
    topo_file.open(topo.c_str());
  }
  catch(std::string s){
    io::messages.add("opening topology failed", "read_input",
		     io::message::error);
    return -1;
  }
  if (!topo_file.is_open()){
    std::cout << "\n\ncould not open " << topo << "!\n" << std::endl;
    io::messages.add("opening topology failed", "read_input",
		     io::message::error);
    return -1;
  }

  in_topo.stream(topo_file);
  in_topo.read(sim.topo, sim.sim.param());

  if(pttopo != ""){
    
    try{
      pttopo_file.open(pttopo.c_str());
    }
    catch(std::string s){
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::error);
      return -1;
    }
    if (!pttopo_file.is_open()){
      std::cout << "\n\ncould not open " << pttopo << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::error);
      return -1;
    }
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.quiet = quiet;
    ipt.read(sim.topo, sim.sim.param());

  }

  // do this after reading in a perturbation topology
  sim.sim.multibath().calculate_degrees_of_freedom(sim.topo);

  if (conf != ""){
    try{
      conf_file.open(conf.c_str());
    }
    catch(std::string s){
      io::messages.add("opening configuration failed", "read_input",
		       io::message::error);
      return -1;
    }
    if (!conf_file.is_open()){
      std::cout << "\n\ncould not open " << conf << "!\n" << std::endl;
      io::messages.add("opening configuration failed", "read_input",
		       io::message::error);
      return -1;
    }

    io::In_Configuration ic(conf_file);
    ic.quiet = quiet;
    ic.read(sim.conf, sim.topo, sim.sim.param());
    sim.conf.initialise(sim.topo, sim.sim.param());
    
  }

  return 0;
}



