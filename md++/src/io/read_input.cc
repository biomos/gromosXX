/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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

#include <util/coding.h>

#ifdef HAVE_HOOMD
#include <HOOMD_GROMOSXX_processor.h>
#endif
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

  // error if no perturbed parameters were read from pttop or restraints
  if(!sim.param().perturbation.perturbed_par && sim.param().perturbation.perturbation){
      io::messages.add("Neither perturbed restraints nor perturbed topology found - if you do not want to perturb anything, turn off PERTURBATION",
		       "read_input", io::message::error);
      return -1;
  }

  sim.multibath().calculate_degrees_of_freedom(topo,
            sim.param().rottrans.rottrans,
            sim.param().posrest.posrest == simulation::posrest_const,
            sim.param().boundary.dof_to_subtract,
            sim.param().dihrest.dihrest == simulation::dihedral_constr,
            sim.param().angrest.angrest == simulation::angle_constr);

  // check the bath parameters
  sim.multibath().check_state(topo.num_atoms());

  if (read_configuration(args, topo, conf, sim, os, quiet) != 0) return -1;

#ifdef HAVE_HOOMD
  // create HOOMD Processor after input files read in successfully
  switch (sim.param().hoomd.processor) {
    case simulation::cpu: sim.proc = boost::shared_ptr<processor::Processor>(new processor::Processor(processor::CPU)); break;
	case simulation::gpus: sim.proc = boost::shared_ptr<processor::Processor>(new processor::Processor(processor::GPUs)); break;
	default: break;
  }
#endif

  return 0;
}

int io::read_input_repex(io::Argument const & args,
		   topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   algorithm::Algorithm_Sequence & md_seq,
                   int rank,
                   int totalNumberOfThreads,
		   std::ostream & os,  
                   bool quiet)
{
    DEBUG(1, "read_input_repex:Start Parsing:")
    //Read in Parameters
    if (read_parameter(args, sim, os, quiet) != 0){
            std::cerr << "\n\t########################################################\n"
                    << "\n\t\tErrors during initial Parameter reading of IMD!\n"
                    << "\n\t########################################################\n";
            io::messages.display(std::cout);
            io::messages.display(std::cerr);
            std::cout << "ERROR after read_input simID\n";
        return -1;}
    if (check_parameter(sim) != 0){
                    std::cerr << "\n\t########################################################\n"
                    << "\n\t\tErrors during initial checking of Paramters!\n"
                    << "\n\t########################################################\n";
            io::messages.display(std::cout);
            io::messages.display(std::cerr);
            std::cout << "ERROR after check parameter\n";
            return -1;}
 
    //Check MPI Mapping
    DEBUG(1, "read_input_repex:Specific MPI fun:")
    //TODO check if enough threads
    //MPI THREAD SIMULATION SPLITTING.
    unsigned int numReplicas = -1;
    if(sim.param().reeds.reeds == 1){
        numReplicas = sim.param().reeds.num_s;
    }
    else if(sim.param().reeds.reeds == 2){
        numReplicas = sim.param().reeds.num_s *  sim.param().reeds.num_eoff;
    }
    else{
        numReplicas = sim.param().replica.num_l*sim.param().replica.num_T;
    }

    //MPI - MAPPING
    //MPI THREAD SPLITING ONTO Simulation - REPLICAS
    unsigned int threadsPerReplicaSimulation = totalNumberOfThreads / numReplicas;
    unsigned int leftOverThreads = totalNumberOfThreads % numReplicas;
    std::map<unsigned int, unsigned int> thread_id_replica_map; // where is which replica

    unsigned int threadID =0;
    int replica_offset = 0;
    for (unsigned int replicaSimulationID = 0; replicaSimulationID < numReplicas; replicaSimulationID++) {

        for (unsigned int replicaSimulationSubThread = 0; replicaSimulationSubThread < threadsPerReplicaSimulation; replicaSimulationSubThread++) {

            threadID = replicaSimulationSubThread + replicaSimulationID*threadsPerReplicaSimulation+replica_offset;
            thread_id_replica_map.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));

        }
        if(leftOverThreads>0 and replicaSimulationID < leftOverThreads){    //left over threads are evenly distirbuted.

            threadID = threadsPerReplicaSimulation + replicaSimulationID*totalNumberOfThreads+replica_offset;
            thread_id_replica_map.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));
            replica_offset++;
        }    
    }            
    int simulationID = thread_id_replica_map[rank];

    //if enough threads avail
    if (totalNumberOfThreads < numReplicas) {
        if (rank == 0) {
            std::cerr << "\n\t########################################################\n"
                    << "\n\t\tErrors during initial Parameter reading!\n"
                    << "\n\t########################################################\n";
            std::cerr << "\n There were not enough MPI threads assigned to this run!\n"
                    << "FOUND THREADS: " << totalNumberOfThreads << "\tNEED: " << numReplicas << "\n";
            std::cout << "\n There were not enough MPI thread assigned to this run!\n"
                    << "FOUND THREADS: " << totalNumberOfThreads << "\tNEED: " << numReplicas << "\n";
            
        #ifdef XXMPI
            MPI_Finalize();
        #endif
        }
        return -1;
    }

    //initialize for RE-EDS ID dependent parameters.
    if(sim.param().reeds.reeds > 0){
      sim.param().eds=sim.param().reeds.eds_para[simulationID];//choose correct eds informations which are ID dependent. That's why this cannot be done earlier.jaaa?
    }
  
  
    //if any replica Ex block - present   
    if (sim.param().reeds.reeds == false && sim.param().replica.retl == false) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\tErrors during initial Parameter reading! "
                << "\n\t\t    No repex block was satisfied!\n"
                << "\n\t########################################################\n";
        std::cerr << "\n Please add one RE-block (e.g.:REPLICA or REEDS) to the imd file.\n";
        std::cout << "\n Please add one RE-block (e.g.:REPLICA or REEDS)  to the imd file.\n";
        return -1;
    }

  // Read in Topology
  if(!quiet){
    std::cout << std::internal << "\tReading Topology\n";
  }

  if (read_topology(args, topo, sim, md_seq, os, quiet) != 0) return -1;

  if(!quiet){
      std::cout << std::internal << "\tReading Special\n";
  }
  // read this before configuration, as it contains topological data...
  if (read_special(args, topo, conf, sim, os, quiet) != 0) return -1;

  // error if no perturbed parameters were read from pttop or restraints
  if(!sim.param().perturbation.perturbed_par && sim.param().perturbation.perturbation){
      io::messages.add("Neither perturbed restraints nor perturbed topology found - if you do not want to perturb anything, turn off PERTURBATION",
		       "read_input", io::message::error);
      return -1;
  }

  sim.multibath().calculate_degrees_of_freedom(topo,
            sim.param().rottrans.rottrans,
            sim.param().posrest.posrest == simulation::posrest_const,
            sim.param().boundary.dof_to_subtract,
            sim.param().dihrest.dihrest == simulation::dihedral_constr,
            sim.param().angrest.angrest == simulation::angle_constr);

  // check the bath parameters
  sim.multibath().check_state(topo.num_atoms());
  if(!quiet){
    std::cout << std::internal << "\tReading Configuration\n";
    std::cout.flush();
  }

	switch(sim.param().reeds.reeds) {
				case 0:
						numReplicas = sim.param().replica.num_T * sim.param().reeds.num_s;
				case 1:
						numReplicas = sim.param().replica.num_T * sim.param().reeds.num_s;
						break;
				case 2:
						numReplicas = sim.param().reeds.num_eoff * sim.param().reeds.num_s;
						break;
	}

	//check if all coordinate files are present:
	int cont = sim.param().replica.cont;
	if(cont == 1 && rank == 0){
		DEBUG(4, "reading configurations for continous");
		for(int x=0; x<numReplicas; x++ ){

				io::Argument tmpArgs(args);   //copy modified args
				std::multimap< std::string, std::string >::iterator it = tmpArgs.lower_bound(("conf"));
				size_t pos = (*it).second.find_last_of(".");
				std::stringstream tmp;
				tmp << "_" << (x+1);
				(*it).second.insert(pos, tmp.str());

				if(read_configuration(tmpArgs, topo, conf, sim, os, quiet)){
						io::messages.add("\nCould not find coordinate file: "+ std::string(it->second)+"\n\n", io::message::error);
						return -1;
				}
			}
		}
	 else if(rank == 0){
				DEBUG(4, "reading configuration no continous");
				//std::cout << "TEST  " << args.lower_bound(("conf"))->second << "\n";
				if (read_configuration(args, topo, conf, sim, os, quiet) != 0) {
						io::messages.add("Could not find coordinate file: "+ std::string(args.lower_bound(("conf"))->second), io::message::error);
						return -1;
				}
		}

#ifdef HAVE_HOOMD
  // create HOOMD Processor after input files read in successfully
  switch (sim.param().hoomd.processor) {
    case simulation::cpu: sim.proc = boost::shared_ptr<processor::Processor>(new processor::Processor(processor::CPU)); break;
	case simulation::gpus: sim.proc = boost::shared_ptr<processor::Processor>(new processor::Processor(processor::GPUs)); break;
	default: break;
  }
#endif
  if(!quiet){
      io::messages.display(os);
  }

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

  io::In_Parameter ip(input_file);
  ip.quiet = quiet;

  ip.read(sim.param(), os);

  io::messages.add("parameter read from " + args[argname_input] +
          "\n" + util::frame_text(ip.title),
          "read input", io::message::notice);

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
    if (args["print"] == "timings_report")
      sim.param().print.timings_report = true;
  }
   
  // check for errors and abort if there are some
  if (io::messages.contains(io::message::error) ||
      io::messages.contains(io::message::critical))
    return -1;

  // check for replicaExchange
  if (sim.param().replica.retl || sim.param().reeds.reeds) {
      // Check output files
        if( args.count("repout") < 1 )
        {
            io::messages.add("No output file for replica exchange specified! Please provide @repout.\n",
           "read_input", io::message::critical);
          return -1;
        }
        if( args.count("repdat") < 1 )
        {
            io::messages.add("No data file for replica exchange specified! Please provide @repdat.\n",
           "read_input", io::message::critical);
          return -1;
        }

        //Check if any REPEX Block was entered and only one!
        //only one replica Ex block - present
        if(sim.param().reeds.reeds == true && sim.param().replica.retl  == true){
            io::messages.add("\n Please provide only one RE-block in the imd file.\n", "read_input", io::message::critical);
            return -1;
        }
    }
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

  io::In_Topology it(topo_file);
  it.quiet = quiet;

  it.read(topo, sim.param(), os);

  io::messages.add("topology read from " + args[argname_topo] + "\n" + util::frame_text(it.title),
          "read input", io::message::notice);

  // check for errors and about before initialization
  if(io::messages.contains(io::message::error) ||
     io::messages.contains(io::message::critical))
    return -1;


   if(args.count(argname_pttopo)<1 && sim.param().reeds.reeds > 0){
      io::messages.add("REEDS on but no perturbation topology specified",
		       "read_input", io::message::critical);
      return -1;
  }
  else if(args.count(argname_pttopo)<1 && sim.param().eds.eds){
      io::messages.add("EDS on but no perturbation topology specified",
		       "read_input", io::message::critical);
      return -1;
  }

  if(sim.param().perturbation.perturbation || sim.param().eds.eds){
    // if there is no perturbation topology there might still be perturbed
    // distance or df restraints, so only warn and do not abort here --MP
    if(args.count(argname_pttopo)<1){
      io::messages.add("No perturbation topology specified",
		       "read_input", io::message::warning);
    }
    else {

    pttopo_file.open(args[argname_pttopo].c_str());

    if (!pttopo_file.is_open()){
      os << "\n\ncould not open " << args[argname_pttopo] << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::critical);
      return -1;
    }

    io::In_Perturbation ipt(pttopo_file);
    ipt.quiet = quiet;

    ipt.read(topo, sim.param(), os);

    sim.param().perturbation.perturbed_par=true;

    io::messages.add("perturbation topology read from " + args[argname_pttopo] + "\n" + util::frame_text(ipt.title),
		     "read input", io::message::notice);
    }
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

  io::In_Configuration ic(conf_file);
  ic.quiet = quiet;

  ic.read(conf, topo, sim, os);

  io::messages.add("configuration read from " + args[argname_conf] + "\n" + util::frame_text(ic.title),
		   "read input", io::message::notice);

  conf.init(topo, sim.param());

  // check for errors and abort
  if (io::messages.contains(io::message::error) ||
      io::messages.contains(io::message::critical))
    return -1;

  return 0;
}
