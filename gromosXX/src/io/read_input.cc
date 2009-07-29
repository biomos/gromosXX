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
#include <io/parameter/check_parameter.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/create_md_sequence.h>

#include <interaction/forcefield/forcefield.h>

#include <util/replica_data.h>

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
		   std::ostream & os,
		   bool quiet)
{

  if (read_parameter(args, sim, os, quiet) != 0) return -1;
  if (check_parameter(sim) != 0) return -1;

  if (read_topology(args, topo, sim, md_seq, os, quiet) != 0) return -1;
  
  // read this before configuration, as it contains topological data...
  if (read_special(args, topo, conf, sim, os, quiet) != 0) return -1;

  sim.multibath().calculate_degrees_of_freedom
          (topo, sim.param().rottrans.rottrans, sim.param().posrest.posrest == simulation::posrest_const, sim.param().boundary.dof_to_subtract);

  // check the bath parameters
  sim.multibath().check_state(topo.num_atoms());

  if (read_configuration(args, topo, conf, sim, os, quiet) != 0) return -1;
  
  return 0;
}

int io::read_replica_input
(
 io::Argument const & args,
 topology::Topology & topo,
 std::vector<configuration::Configuration> & conf,
 simulation::Simulation & sim,
 algorithm::Algorithm_Sequence & md_seq,
 std::vector<util::Replica_Data> & replica_data,
 std::ostream & os,
 bool quiet)
{
  if (read_parameter(args, sim, os, quiet) != 0) return -1;

  if (read_topology(args, topo, sim, md_seq, os, quiet) != 0) return -1;

  for (unsigned int i = 0; i < conf.size(); ++i) {
    if (read_special(args, topo, conf[i], sim, os, quiet) != 0) return -1;
  }

  sim.multibath().calculate_degrees_of_freedom
          (topo, sim.param().rottrans.rottrans, sim.param().posrest.posrest == simulation::posrest_const, sim.param().boundary.dof_to_subtract);

  // check the bath parameters
  sim.multibath().check_state(topo.num_atoms());

  if (read_replica_configuration(args, topo, conf, sim, replica_data, os, quiet) != 0)
    return -1;


  return 0;
}


int io::read_parameter(io::Argument const & args,
		       simulation::Simulation & sim,
		       std::ostream & os,
		       bool quiet)
{
  io::igzstream input_file;
  
  input_file.open(args[argname_input].c_str());
  
  if (!input_file.is_open()){
    os << "\n\ncould not open " << args[argname_input] << "!\n" << std::endl;
    io::messages.add("opening input failed", "read_input",
		     io::message::critical);
    return -1;
  }
  
  io::messages.add("parameter read from " + args[argname_input],
		   "read input",
		   io::message::notice);
  
  io::In_Parameter ip(input_file);
  ip.quiet = quiet;
  
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
  
  if (args.count("print") > 0){
    if (args["print"] == "pairlist")
      sim.param().pairlist.print = true;
  }

  // check for errors and abort if there are some
  if (io::messages.contains(io::message::error) ||
      io::messages.contains(io::message::critical))
    return -1;
  
  return 0;
}

int io::read_topology(io::Argument const & args,
		      topology::Topology & topo,
		      simulation::Simulation & sim,
		      algorithm::Algorithm_Sequence & md_seq,
		      std::ostream & os,
		      bool quiet)
{
  io::igzstream topo_file, pttopo_file;
  
  topo_file.open(args[argname_topo].c_str());
    if (!topo_file.is_open()){
    os << "\n\ncould not open " << args[argname_topo] << "!\n" << std::endl;
    io::messages.add("opening topology failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("topology read from " + args[argname_topo],
		   "read input",
		   io::message::notice);
  
  io::In_Topology it(topo_file);
  it.quiet = quiet;
  
  it.read(topo, sim.param(), os);

  // check for errors and about before initialization
  if(io::messages.contains(io::message::error) ||
     io::messages.contains(io::message::critical))
    return -1;
  
  if(sim.param().perturbation.perturbation || sim.param().eds.eds){
    if(args.count(argname_pttopo)<1){
      io::messages.add("No perturbation topology specified",
		       "read_input", io::message::critical);
      return -1;
    }
    
    pttopo_file.open(args[argname_pttopo].c_str());
    
    if (!pttopo_file.is_open()){
      os << "\n\ncould not open " << args[argname_pttopo] << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::critical);
      return -1;
    }
    
    io::messages.add("perturbation topology read from " + args[argname_pttopo],
		     "read input",
		     io::message::notice);
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.quiet = quiet;
    
    ipt.read(topo, sim.param());
    
  }
  
  topo.init(sim, os, quiet);

  // and create the algorithms
  // (among them the forcefield!)
  algorithm::create_md_sequence(md_seq, topo, sim, it, os, quiet);

  return 0;
}

int io::read_configuration(io::Argument const & args,
			   topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   std::ostream & os,
			   bool quiet)
{
  io::igzstream conf_file;

  DEBUG(7, "reading configuration");
  conf_file.open(args[argname_conf].c_str());

  if (!conf_file.is_open()){
    os << "\n\ncould not open " << args[argname_conf] << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("configuration read from " + args[argname_conf],
		   "read input",
		   io::message::notice);
  
  io::In_Configuration ic(conf_file);
  ic.quiet = quiet;
  
  ic.read(conf, topo, sim, os);

  conf.init(topo, sim.param());

  // check for errors and abort
  if (io::messages.contains(io::message::error) || 
      io::messages.contains(io::message::critical))
    return -1;
    
  return 0;
}

/**
 * read in a configuration
 */
int io::read_replica_configuration
(
 io::Argument const &args,
 topology::Topology &topo,
 std::vector<configuration::Configuration> & conf,
 simulation::Simulation & sim,
 std::vector<util::Replica_Data> & replica_data,
 std::ostream & os,
 bool quiet
 )
{
  io::igzstream conf_file;

  DEBUG(7, "reading replica configurations");
  conf_file.open(args[argname_conf].c_str());

  if (!conf_file.is_open()){
    os << "\n\ncould not open " << args[argname_conf] << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
    return -1;
  }

  io::messages.add("replica configurations read from " + args[argname_conf],
		   "read input", io::message::notice);
  
  io::In_Configuration ic(conf_file);
  ic.quiet = quiet;
  
  ic.read_replica(conf, topo, sim, replica_data, os);

  for(unsigned int i=0; i<conf.size(); ++i)
    conf[i].init(topo, sim.param());

  // check for errors and abort
  if (io::messages.contains(io::message::error) || 
      io::messages.contains(io::message::critical))
    return -1;
    
  return 0;
}
