/**
 * @file vtest.cc
 * simple simulation program.
 */

#include <iostream>
#include <iomanip>

#include <math/gmath.h>
#include <simulation/simulation.h>
#include <interaction/interaction.h>
#include <algorithm/algorithm.h>

// sppecial includes
#include <algorithm/integration/runge_kutta.h>

using namespace math;

int main()
{
  const int SIZE = 2;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.num_solute_atoms(SIZE);
  
  // initialize
  the_system.pos()(0) = Vec( -1.0, 0.0, 0.0 );
  the_system.pos()(1) = Vec(  1.0, 0.0, 0.0 );
  
  the_system.vel() = 0.0;
  the_system.force() = 0.0;
  
  the_topology.mass() = 1.0;
  
  // simulation
  typedef simulation::simulation<simulation::topology,
    simulation::system> simulation_type;
  
  simulation_type the_simulation(the_topology, the_system);
  
  // forcefield
  interaction::forcefield<simulation_type> the_forcefield;
  
  interaction::harmonic_bond_interaction<simulation_type> 
    *the_bond_interaction =
    new interaction::harmonic_bond_interaction<simulation_type>;
  
  // add a bond type and a bond
  the_bond_interaction->add(1000, 1.5);
  the_simulation.topology().bonds().add(0, 1, 0);
  
  // add to the forcefield
  the_forcefield.add_interaction(the_bond_interaction);
  
  // create the algorithm
  algorithm::runge_kutta<simulation_type> RK;

  // simulate
  for(int i=0; i<10000; ++i){
    // std::cout << "STEP " << i << std::endl;
    // std::cout << "POSITION" << std::endl;
    // std::cout << the_simulation.system().pos() << std::endl;
    // std::cout << "VELOCITY" << std::endl;
    // std::cout << the_simulation.system().vel() << std::endl;

    std::cout << std::setw(8) << i << std::setw(20)
	      << the_simulation.system().pos()(1)(0) 
	      << std::setw(20) << the_simulation.system().vel()(1)(0)
	      << std::endl;
    
    RK.step(the_simulation, the_forcefield, 0.001);
    // algorithm::leap_frog<simulation_type>
    // ::step(the_simulation, the_forcefield, 0.001);
    
  }
  
  std::cout << "VTEST finished successfully" << std::endl;
  
  return 0;
}

    
