/**
 * @file vtest.cc
 * simple simulation program.
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

// helper functions to prepare the md run
#include "create_forcefield.h"

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

    // parse print flag
    int print_pairlist = 0, print_force = 1;
    { // print
      io::Argument::const_iterator it = args.lower_bound("print"),
	to = args.upper_bound("print");
      if (it != to)
	std::cout << "printing\n";
      for( ; it != to; ++it){
	std::string s;
	int num;
	std::string::size_type sep = it->second.find(':');
	if (sep == std::string::npos){
	  s = it->second;
	  num = 1;
	}
	else{
	  s = it->second.substr(0, sep);
	  num = atoi(it->second.substr(sep+1, std::string::npos).c_str());
	}
	std::cout << "\t" << std::setw(15) << s << std::setw(6) << num << "\n";

	if (s == "pairlist") print_pairlist = num;
	else if (s == "force") print_force = num;
	else throw std::string("unknown @print argument");
      }
    }

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
  
    simulation::System<math::any> the_system;
    simulation::Topology the_topology;

    DEBUG(7, "reading the files");

    // read in the files - those are necessary
    std::ifstream topo_file(args["topo"].c_str());
    io::InTopology topo(topo_file);
  
    std::ifstream sys_file(args["struct"].c_str());
    io::InTrajectory sys(sys_file);

    std::ifstream input_file(args["input"].c_str());
    io::InInput input(input_file);

    topo >> the_topology;
    sys >> the_system;

    DEBUG(7, "topology and system read");
  
    // simulation
    typedef simulation::Simulation<simulation::Topology,
      simulation::System<math::any> > simulation_type;
  
    simulation_type the_simulation(the_topology, the_system);

    typedef interaction::twin_range_pairlist_cg<simulation_type> pairlist_type;

    // add solvent
    int nsm;
    input.read_SYSTEM(nsm);
    if (nsm) the_simulation.solvate(0, nsm);

    // FORCEFIELD
    interaction::Forcefield<simulation_type> &the_forcefield
      = create_forcefield<simulation_type, pairlist_type>(topo, input, the_topology);

    input >> the_simulation;
  
    // create the algorithm
    int ntc;
    double tolerance;
    input.read_SHAKE(ntc, tolerance);

    algorithm::runge_kutta<simulation_type> RK;
    algorithm::Leap_Frog<simulation_type> LF;
    algorithm::Shake<simulation_type> shake(tolerance);

    // prepare for the run
    the_simulation.calculate_degrees_of_freedom();
    
    // prepare for the output
    int print_trajectory, print_velocity, print_energy;
    input.read_PRINT(print_trajectory, print_velocity, print_energy);

    std::ofstream trap(args["trj"].c_str());  // trajectory is required
    std::ofstream final(args["fin"].c_str()); // final structure is required

    // use a G96 trajectory
    io::OutG96Trajectory<simulation_type> traj(trap, final, print_trajectory);

    // optional files
    std::ofstream trav; // velocity trajectory
    if (args.count("trv") == 1){
      trav.open(args["trv"].c_str());
      traj.velocity_trajectory(trav, print_velocity);
    }
    
    std::ofstream traf; // force trajectory
    if (args.count("trf") == 1){
      traf.open(args["trf"].c_str());
      traj.force_trajectory(traf, print_force);
    }

    std::ostream *trprint = &std::cout;
    if (args.count("trp") == 1){
      trprint = new std::ofstream(args["trp"].c_str());
    }
    
    traj.print_title("\tvtest(gromosXX) MD simulation");

    std::cout << the_simulation.multibath();

    std::cout << "Messages (startup)\n";
    if (io::messages.display(std::cout) > io::message::warning)
      return 1;
    std::cout << "\n";
    io::messages.clear();

    int num_steps;
    double t0, dt;
    input.read_STEP(num_steps, t0, dt);
    std::cout << "steps: " << num_steps << " dt: " << dt << std::endl;

    the_simulation.time(t0);

    std::cout << "starting MD\n\n";

    // simulate
    for(int i=0; i<num_steps; ++i){
    
      traj << the_simulation;

      if (runge_kutta)
	RK.step(the_simulation, the_forcefield, dt);
      else // leap frog
	LF.step(the_simulation, the_forcefield, dt);
  
      if (print_energy && the_simulation.steps() % print_energy == 0){
	std::cout << the_simulation.multibath();
      }

      if (print_pairlist && the_simulation.steps() % print_pairlist == 0){

	std::vector<interaction::Interaction<simulation_type> *>::const_iterator it = the_forcefield.begin(),
	  to = the_forcefield.end();
	
	for( ; it != to; ++it){
	  
	  if ((*it)->name == "NonBonded"){

	    *trprint << "shortrange\n" 
		     << dynamic_cast<interaction::Nonbonded_Interaction<simulation_type,
	      pairlist_type> *>(*it)->pairlist().short_range()
		     << std::endl;
	    *trprint << "longrange\n" 
		     << dynamic_cast<interaction::Nonbonded_Interaction<simulation_type,
	      pairlist_type> *>(*it)->pairlist().long_range()
		     << std::endl;
	  }
	  
	}
	
      }
      
      try{
	std::cout << "shake solute:  " << shake.solute(the_topology, the_system, dt)
		  << "\n";
	std::cout << "shake solvent: " << shake.solvent(the_topology, the_system, dt)
		  << "\n";
      }
      catch(std::runtime_error e){
	the_system.exchange_pos();
	traj << the_simulation;
	throw;
      }
    
      the_simulation.increase_time(dt);
    
    }

    std::cout << "\nwriting final structure" << std::endl;

    traj << io::final << the_simulation;
  
    std::cout << "\nVTEST finished successfully\n\n" << std::endl;
  
    std::cout << "messages (simulation)\n";
    io::messages.display(std::cout);
    std::cout << "\n\n";
  
  }
  catch(std::runtime_error e){
    // runtime problem
    std::cout << "severe error encountered:\n"
	      << e.what() << std::endl;
    
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

