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
    
  // check command line
  if (argc < 4){
    std::cout << "usage: " << argv[0] << " topology structure input "
      "[RungeKutta] [debuglevel]\n\n";
    io::messages.add("wrong number of arguments.", 
		     "vtest", 
		     io::message::critical);
    
    throw std::runtime_error("wrong arguments");
  }
  
  bool runge_kutta = false;
  if (argc > 4){
    std::string s = argv[4];
    if (s == "RungeKutta"){
      runge_kutta = true;
      io::messages.add("using Runge Kutta integration scheme",
		       "vtest",io::message::notice);
    }
    else{
      io::messages.add("using Leap Frog integration scheme",
		       "vtest", io::message::notice);
    }
  }
  if (argc > 5){
    debug_level = atoi(argv[5]);
    std::cout << "setting debug_level = " << debug_level << std::endl;
  }
  
  simulation::system the_system;
  simulation::Topology the_topology;

  DEBUG(7, "reading the files");
  
  // read in the files
  std::ifstream topo_file(argv[1]);
  io::InTopology topo(topo_file);
  
  std::ifstream sys_file(argv[2]);
  io::InTrajectory sys(sys_file);

  std::ifstream input_file(argv[3]);
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
  
  // nonbonded
  typedef interaction::twin_range_pairlist_cg<simulation_type> pairlist_type;

  interaction::Nonbonded_Interaction<simulation_type, pairlist_type>
    *the_nonbonded_interaction =
    new interaction::Nonbonded_Interaction<simulation_type, pairlist_type>;

  DEBUG(7, "parsing parameter");

  // read parameter
  topo >> *the_bond_interaction;
  topo >> *the_angle_interaction;
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
  std::ofstream trap("vtest.trj");
  std::ofstream trav("vtest.trv");
  std::ofstream traf("vtest.trf");
  std::ofstream final("vtest.fin");
  
  io::OutTrajectory<simulation_type> traj(trap, final);
  traj.velocity_trajectory(trav);
  traj.force_trajectory(traf);

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
  
    std::cout << "shortrange\n" 
	      << the_nonbonded_interaction->pairlist().short_range()
	      << std::endl;
    std::cout << "longrange\n" 
	      << the_nonbonded_interaction->pairlist().long_range()
	      << std::endl;

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
    std::cout << "severe error encountered:\n";
    io::messages.display();
    std::cout << std::endl;
    
    return 1;
  }
  
  return 0;
}

