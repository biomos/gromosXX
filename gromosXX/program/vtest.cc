/**
 * @file vtest.cc
 * simple simulation program.
 */

#include <iostream>
#include <iomanip>
#include <fstream>

#include <math/gmath.h>
#include <simulation/simulation.h>
#include <interaction/interaction.h>
#include <algorithm/algorithm.h>
#include <io/io.h>

// sppecial includes
#include <algorithm/integration/runge_kutta.h>

using namespace math;

int main(int argc, char *argv[])
{
  try{
    
  // check command line
  if (argc != 3){
    std::cout << "usage: " << argv[0] << " topology structure\n\n";
    throw std::runtime_error("wrong arguments");
  }
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  // read in the files
  std::ifstream topo_file(argv[1]);
  io::InTopology topo(topo_file);
  
  std::ifstream sys_file(argv[2]);
  io::InTrajectory sys(sys_file);

  topo >> the_topology;
  sys >> the_system;
  
  // simulation
  typedef simulation::simulation<simulation::topology,
    simulation::system> simulation_type;
  
  simulation_type the_simulation(the_topology, the_system);
  
  // forcefield
  interaction::forcefield<simulation_type> the_forcefield;
  
  interaction::harmonic_bond_interaction<simulation_type> 
    *the_bond_interaction =
    new interaction::harmonic_bond_interaction<simulation_type>;
  
  // read parameter
  topo >> *the_bond_interaction;
  
  // add to the forcefield
  the_forcefield.add_interaction(the_bond_interaction);
  
  // create the algorithm
  algorithm::runge_kutta<simulation_type> RK;

  // prepare for the output
  std::ofstream trap("vtest.trj");
  std::ofstream trav("vtest.trv");
  std::ofstream final("vtest.fin");
  
  io::OutTrajectory<simulation_type> traj(trap, final);
  traj.velocity_trajectory(trav);

  if (io::messages.display() > io::message::notice)
    return 1;
  
  // simulate
  for(int i=0; i<100; ++i){
    
    traj << the_simulation;

    RK.step(the_simulation, the_forcefield, 0.001);
    // algorithm::leap_frog<simulation_type>
    // ::step(the_simulation, the_forcefield, 0.001);
    
  }

  traj << io::final << the_simulation;
  
  std::cout << "VTEST finished successfully" << std::endl;
  
  io::messages.display();
  
  }
  catch(std::runtime_error e){
    io::messages.display();
    return 1;
  }
  
  return 0;
}

