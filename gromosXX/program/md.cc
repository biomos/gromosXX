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
// #include <io/trajectory/InFlexibleConstraints.h>

// sppecial includes
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
	"trj", "fin", "trv", "trf", "tre", "print", "trp"
      };
    
    int nknowns = 13;
    
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

    // determine whether we do perturbation
    bool perturbation = false;
    if (args.count("pert") == 1){
      perturbation = true;
      if (runge_kutta)
	io::messages.add("perturbation with runge kutta integration"
			 " not allowed",
			 "md", io::message::error);
    }
    
    // this is not the nicest solution...
    if (runge_kutta){

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
      
      if (the_MD.initialize(args)){
	return 1;
      }
      
      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

    }
    else if (perturbation){ // leap frog + perturbation
      // topology and system
      simulation::System<math::any> the_system;
      simulation::Perturbation_Topology the_topology;

      // simulation
      typedef simulation::Simulation<simulation::Perturbation_Topology,
	simulation::System<math::any> > simulation_type;
  
      simulation_type the_simulation(the_topology, the_system);
      
      algorithm::Perturbation_MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Shake<simulation_type>,
	algorithm::Leap_Frog<simulation_type> >
	the_MD(the_simulation);
    
      if(the_MD.initialize(args)){
	return 1;
      }

      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

    }
    else{ // leap frog, no perturbation
      // topology and system
      simulation::System<math::any> the_system;
      simulation::Topology the_topology;

      // simulation
      typedef simulation::Simulation<simulation::Topology,
	simulation::System<math::any> > simulation_type;
  
      simulation_type the_simulation(the_topology, the_system);
      
      algorithm::MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Shake<simulation_type>,
	algorithm::Leap_Frog<simulation_type> >
	the_MD(the_simulation);
    
      if(the_MD.initialize(args)){
	return 1;
      }

      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

      simulation::Energy energy, fluctuation;
      the_MD.simulation().system().energy_averages().average(energy, fluctuation);
      
      io::print_ENERGY(std::cout, energy,
		       the_MD.simulation().topology().energy_groups(), "AVERAGE ENERGIES");

      io::print_ENERGY(std::cout, fluctuation,
		       the_MD.simulation().topology().energy_groups(), "ENERGY FLUCTUATIONS");

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

