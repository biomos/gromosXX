/**
 * @file md.cc
 * the md program.
 */

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include "../src/debug.h"

#include <math/gmath.h>
#include <io/message.h>
#include <simulation/core.h>
#include <math/periodicity.h>
#include <simulation/simulation.h>
#include <simulation/perturbation.h>
#include <interaction/interaction.h>
#include <io/io.h>
#include <algorithm/algorithm.h>

// special includes
#include <algorithm/integration/runge_kutta.h>

// global variables for debug
#include "../src/debug.cc"

using namespace math;

int main(int argc, char *argv[])
{
  try{
    
    char *knowns[] = 
      {
	"topo", "struct", "input", "verb", "alg", "pert",
	"trj", "fin", "trv", "trf", "tre", "trg", "print", "trp"
      };
    
    int nknowns = 14;
    
    string usage = argv[0];
    usage += "\n\t@topo    <topology>\n";
    usage += "\t[@pert   <perturbation topology>]\n";
    usage += "\t@struct  <coordinates>\n";
    usage += "\t@input   <input>\n";
    usage += "\t@trj     <trajectory>\n";
    usage += "\t@fin     <final structure>\n";
    usage += "\t[@trv    <velocity trajectory>]\n";
    usage += "\t[@trf    <force trajectory>]\n";
    usage += "\t[@tre    <energy trajectory>]\n";
    usage += "\t[@trg    <free energy trajectory>]\n";
    usage += "\t[@alg    <RK|LF>]\n";
    usage += "\t[@print  <pairlist/force>]\n";
    usage += "\t[@trp    <print file>]\n";
    usage += "\t[@verb   <[module:][submodule:]level>]\n";

    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    parse_verbosity(args);

    // determine which algorithm to use
    bool runge_kutta = false;
    if (args.count("alg") != -1){
      if (args["alg"] == "RK"){
	runge_kutta = true;
	io::messages.add("using Runge Kutta integration scheme",
			 "md",io::message::notice);
      }
      else if(args["alg"] == "LF"){
	io::messages.add("using Leap Frog integration scheme",
			 "md",io::message::notice);
      }
      else{
	io::messages.add("unknown integration scheme (@alg) " + args["alg"],
			 "md",io::message::error);
      }
    }
    
    // determine whether we do perturbation and RK
    bool perturbation = false;
    if (args.count("pert") == 1)
      perturbation = true;

    if (runge_kutta){

      if (perturbation){
	io::messages.add("perturbation with runge kutta integration"
			 " not allowed",
			 "md", io::message::error);
      }

      /**
      // topology and system
      simulation::System<math::any> the_system;
      simulation::Perturbation_Topology the_topology;

      // simulation
      typedef simulation::Simulation<simulation::Perturbation_Topology,
	simulation::System<math::any> > simulation_type;
  
      simulation_type the_simulation(the_topology, the_system);

      algorithm::MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Shake<simulation_type>,
	algorithm::runge_kutta<simulation_type> >
	the_MD(the_simulation);
      
      if (do_md(the_MD, args)){
	return 1;
      }      
      */
      return 1;
    }
    else if (perturbation){ // leap frog + perturbation

      algorithm::Perturbation_MD<
	algorithm::perturbed_MD_spec
	> 
	the_MD;

      if (the_MD.do_md(args)){
	return 1;
      }
    }
    else{ // leap frog, no perturbation
      algorithm::MD<
	algorithm::MD_spec
	> 
	the_MD;

      if (the_MD.do_md(args)){
	return 1;
      }
    }
    
    std::cout << "\nMD finished successfully\n\n" << std::endl;
  
    std::cout << "messages (simulation)\n";
    io::messages.display(std::cout);
    std::cout << "\n\n";
  
  }
  catch(std::runtime_error e){
    // runtime problem
    std::cout << "severe error encountered:\n"
	      << e.what() << std::endl;
    std::cout << "messages (simulation)\n";
    
    io::messages.display();
    std::cout << std::endl;
    
    return 1;
  }
  catch(std::string s){
    // argument exception
    std::cout << s << std::endl;
    return 2;
  }
  
  return 0;
}

