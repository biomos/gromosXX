/**
 * @file blockinput.t.cc
 * test routines for blockinput.
 */

#include "../math/gmath.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <stdexcept>

#include "../simulation/simulation.h"
#include "../interaction/interaction.h"

#include "io.h"

/**
 * test the io.
 */
int test_blockio()
{
  int result = 0;

  simulation::system the_system;
  simulation::topology the_topology;
  
  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  simulation_type the_simulation(the_topology, the_system);

  interaction::forcefield<simulation_type> the_forcefield;

  interaction::harmonic_bond_interaction<simulation_type> *bond_interaction
    = new interaction::harmonic_bond_interaction<simulation_type>;
  
  std::ifstream topo_file("/home/markus/test/hexa10.topo");
  io::InTopology topo(topo_file);
  
  std::ifstream sys_file("/home/markus/test/hexa10.coord");
  io::InTrajectory sys(sys_file);

  topo >> *bond_interaction;
  topo >> the_topology;

  sys >> the_system;

  if (io::message::notice != io::messages.display()) ++result;

  // output
  std::ofstream final("blockinput.t.fin");
  std::ofstream trj("blockinput.t.trj");
  
  io::OutTrajectory<simulation_type> traj(trj, final);
  traj << the_simulation;
  // traj << io::decorated << the_simulation;
  traj << io::final << the_simulation;

  return result;

}

int main()
{
  int r1;
  if ((r1 = test_blockio()))
    std::cout << "test_blockio failed" << std::endl;
  
  return r1;
}

