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

    // parse print flag
    int print_pairlist = 0, print_force = 1;
    { // print
      io::Argument::const_iterator it = args.lower_bound("print"),
	to = args.upper_bound("print");
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
  
    simulation::system the_system;
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

    int nsm;
    input.read_SYSTEM(nsm);
    if (nsm) the_topology.solvate(0, nsm);
  
    // simulation
    typedef simulation::Simulation<simulation::Topology,
      simulation::system> simulation_type;
  
    simulation_type the_simulation(the_topology, the_system);
  
    // FORCEFIELD
    interaction::forcefield<simulation_type> the_forcefield;

    // bonds
    interaction::harmonic_bond_interaction<simulation_type> 
      *the_bond_interaction =
      new interaction::harmonic_bond_interaction<simulation_type>;

    // angles
    interaction::angle_interaction<simulation_type>
      *the_angle_interaction = 
      new interaction::angle_interaction<simulation_type>;
  
    // improper dihedrals
    interaction::Improper_dihedral_interaction<simulation_type>
      *the_improper_interaction = 
      new interaction::Improper_dihedral_interaction<simulation_type>;
    
    // dihedrals
    interaction::Dihedral_interaction<simulation_type>
      *the_dihedral_interaction =
      new interaction::Dihedral_interaction<simulation_type>;
    
    // nonbonded
    typedef interaction::twin_range_pairlist_cg<simulation_type> pairlist_type;

    interaction::Nonbonded_Interaction<simulation_type, pairlist_type>
      *the_nonbonded_interaction =
      new interaction::Nonbonded_Interaction<simulation_type, pairlist_type>;

    DEBUG(7, "parsing parameter");

    // read parameter
    topo >> *the_bond_interaction;
    topo >> *the_angle_interaction;
    topo >> *the_improper_interaction;
    topo >> *the_nonbonded_interaction;
  
    input >> the_simulation;

    // add to the forcefield
    bool do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
    input.read_FORCE(do_bond, do_angle, do_improper,
		     do_dihedral, do_nonbonded);
  
    if (do_bond)
      the_forcefield.add_interaction(the_bond_interaction);
    if (do_angle)
      the_forcefield.add_interaction(the_angle_interaction);
    if (do_improper)
      the_forcefield.add_interaction(the_improper_interaction);
    if (do_dihedral)
      the_forcefield.add_interaction(the_dihedral_interaction);
    if (do_nonbonded)
      the_forcefield.add_interaction(the_nonbonded_interaction);

    // decide on SHAKE
    int ntc;
    double tolerance;
    input.read_SHAKE(ntc, tolerance);

    switch(ntc){
      case 1:
	break;
      case 2: 
	std::cout << "SHAKE bonds containing hydrogen atoms" << std::endl;
	the_topology.solute().
	  add_bond_length_constraints(1.0,
				      the_topology.mass(),
				      the_bond_interaction->parameter());
	break;
      case 3: 
	std::cout << "SHAKE all bonds" << std::endl;
	the_topology.solute().
	  add_bond_length_constraints(the_bond_interaction->parameter());
	break;
      default:
	std::cout << "wrong ntc" << std::endl;
    }
  
    // create the algorithm
    algorithm::runge_kutta<simulation_type> RK;
    algorithm::Shake<simulation_type> shake(tolerance);

    // prepare for the output
    int print_trajectory, print_velocity, print_energy;
    input.read_PRINT(print_trajectory, print_velocity, print_energy);

    std::ofstream trap(args["trj"].c_str());  // trajectory is required
    std::ofstream final(args["fin"].c_str()); // final structure is required

    io::OutTrajectory<simulation_type> traj(trap, final, print_trajectory);

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

    std::cout << "Messages (startup)\n";
    if (io::messages.display(std::cout) > io::message::notice)
      return 1;
    std::cout << "\n";

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
      else
	algorithm::leap_frog<simulation_type>
	  ::step(the_simulation, the_forcefield, dt);
  
      if (print_pairlist && the_simulation.steps() % print_pairlist == 0){
	*trprint << "shortrange\n" 
		 << the_nonbonded_interaction->pairlist().short_range()
		 << std::endl;
	*trprint << "longrange\n" 
		 << the_nonbonded_interaction->pairlist().long_range()
		 << std::endl;
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

