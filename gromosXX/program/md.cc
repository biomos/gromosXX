/**
 * @file md.cc
 * the md program.
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include "../src/debug.h"

#include <math/gmath.h>
#include <io/message.h>
#include <simulation/simulation.h>
#include <interaction/interaction.h>
#include <io/io.h>
#include <algorithm/algorithm.h>

// sppecial includes
#include <algorithm/integration/runge_kutta.h>

#include "../src/debug.cc"

using namespace math;

int main(int argc, char *argv[])
{
  try{
    
    char *knowns[] = 
      {
	"topo", "struct", "input", "verb", "alg",
	"trj", "fin", "trv", "trf", "print", "trp"
      };
    
    int nknowns = 11;
    
    string usage = argv[0];
    usage += "\n\t@topo    <topology>\n";
    usage += "\t@struct  <coordinates>\n";
    usage += "\t@input   <input>\n";
    usage += "\t@trj     <trajectory>\n";
    usage += "\t@fin     <final structure>\n";
    usage += "\t@trv     <velocity trajectory>\n";
    usage += "\t@trf     <force trajectory>\n";
    usage += "\t@alg     <RK|LF>\n";
    usage += "\t@print   <pairlist>\n";
    usage += "\t@trp     <print file>\n";
    usage += "\t@verb    <[module:][submodule:]level>\n";
    
    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    parse_verbosity(args);

    // determine which algorithm to use
    bool runge_kutta = false;
    if (args.count("alg") != -1){
      if (args["alg"] == "RK"){
	runge_kutta = true;
	io::messages.add("using Runge Kutta integration scheme",
			 "vtest",io::message::notice);
      }
      else if(args["alg"] == "LF"){
	io::messages.add("using Leap Frog integration scheme",
			 "vtest",io::message::notice);
      }
      else{
	io::messages.add("unknown integration scheme (@alg) " + args["alg"],
			 "vtest",io::message::error);
      }
    }

    
    // topology and system
    simulation::System<math::any> the_system;
    simulation::Topology the_topology;

    // simulation
    typedef simulation::Simulation<simulation::Topology,
      simulation::System<math::any> > simulation_type;
  
    simulation_type the_simulation(the_topology, the_system);

    // this is not the nicest solution...
    if (runge_kutta){

      algorithm::MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Shake<simulation_type>,
	algorithm::runge_kutta<simulation_type> >
	the_MD(the_simulation);
      
      the_MD.initialize(args);
      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

    }
    else{
      
      algorithm::MD<simulation_type,
	algorithm::Berendsen_Thermostat,
	algorithm::Berendsen_Barostat,
	algorithm::Shake<simulation_type>,
	algorithm::Leap_Frog<simulation_type> >
	the_MD(the_simulation);
      
      the_MD.initialize(args);
      the_MD.run();

      std::cout << "\nwriting final structure" << std::endl;
      the_MD.trajectory() << io::final << the_MD.simulation();

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

