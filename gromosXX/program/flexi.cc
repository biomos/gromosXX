/**
 * @file flexi.cc
 * the md program with flexible constraints.
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
	"topo", "struct", "input", "verb", "pert",
	"trj", "fin", "trv", "trf", "tre", "print", "trp",
	"flexcon", "trc"
      };
    
    int nknowns = 14;
    
    string usage = argv[0];
    usage += "\n\t@topo     <topology>\n";
    usage += "\t[@pert    <perturbation topology>]\n";
    usage += "\t@struct   <coordinates>\n";
    usage += "\t@input    <input>\n";
    usage += "\t@trj      <trajectory>\n";
    usage += "\t@fin      <final structure>\n";
    usage += "\t[@trv     <velocity trajectory>]\n";
    usage += "\t[@trf     <force trajectory>]\n";
    usage += "\t[@tre     <energy trajectory>]\n";
    usage += "\t[@print   <pairlist/force>]\n";
    usage += "\t[@trp     <print file>]\n";
    usage += "\t[@flexcon <flexible constraints data file>]\n";
    usage += "\t[@trc     <flexible constraints output file>]\n";
    usage += "\t[@verb    <[module:][submodule:]level>]\n";

    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    parse_verbosity(args);

    // determine whether we do perturbation
    bool perturbation = false;
    if (args.count("pert") == 1){
      perturbation = true;
    }
    
    if (perturbation){ // leap frog + perturbation
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
	algorithm::Perturbed_Flexible_Constraint<simulation_type>,
	algorithm::Leap_Frog<simulation_type> >
	the_MD(the_simulation);
    
      if(the_MD.initialize(args)){
	return 1;
      }

      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

      if (args.count("trc") == 1){
	std::cout << "\nwriting flexible constraints data" << std::endl;
	std::ofstream constr_file(args["trc"].c_str());
	io::OutFlexibleConstraints flexout(constr_file);
	flexout.write_title(the_MD.title);
	flexout.write_FLEXCON(the_MD.distance_constraint_algorithm().vel(),
			      the_MD.simulation().topology());

      }
      else{
	std::cout << "writing of final flexible constraint data"
		  << " not required\n";
	
      }
      
      simulation::Energy energy, fluctuation;
      the_MD.simulation().system().energy_averages().average(energy, fluctuation);
      
      io::print_ENERGY(std::cout, energy,
		       the_MD.simulation().topology().energy_groups(),
		       "AVERAGE ENERGIES");

      io::print_ENERGY(std::cout, fluctuation,
		       the_MD.simulation().topology().energy_groups(),
		       "ENERGY FLUCTUATIONS");


    }
    else{ // leap frog, no perturbation
      // topology and system
      std::cout << "no perturbation" << std::endl;
      simulation::System<math::any> the_system;
      simulation::Topology the_topology;

      std::cout << "system & topology" << std::endl;
      
      // simulation
      typedef simulation::Simulation<simulation::Topology,
	simulation::System<math::any> > simulation_type;
  
      simulation_type the_simulation(the_topology, the_system);
      
      std::cout << "simulation" << std::endl;
      
      algorithm::MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Flexible_Constraint<simulation_type>,
	algorithm::Leap_Frog<simulation_type> >
	the_MD(the_simulation);
    
      std::cout << "the md algorithm" << std::endl;
      
      if(the_MD.initialize(args)){
	return 1;
      }
      std::cout << "initialized" << std::endl;
      
      the_MD.run();
      
      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

      if (args.count("trc") == 1){
	std::cout << "writing flexible constraints data\n\n";
	std::ofstream constr_file(args["trc"].c_str());
	io::OutFlexibleConstraints flexout(constr_file);
	flexout.write_title(the_MD.title);
	flexout.write_FLEXCON(the_MD.distance_constraint_algorithm().vel(),
			      the_MD.simulation().topology());

      }
      else{
	std::cout << "writing of final flexible constraint data"
		  << " not required\n";
	// std::cout << "trc: " << args.count("trc") << std::endl;
	
      }

      simulation::Energy energy, fluctuation;
      the_MD.simulation().system().energy_averages().average(energy, fluctuation);
      
      io::print_ENERGY(std::cout, energy,
		       the_MD.simulation().topology().energy_groups(),
		       "AVERAGE ENERGIES");

      io::print_ENERGY(std::cout, fluctuation,
		       the_MD.simulation().topology().energy_groups(),
		       "ENERGY FLUCTUATIONS");

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

