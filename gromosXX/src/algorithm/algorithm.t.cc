/**
 * @file algorithm.t.cc
 * test routines for the algorithm library
 */

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "../math/gmath.h"

#include "../simulation/simulation.h"

#include "../interaction/interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "integration/leap_frog.h"
#include "integration/runge_kutta.h"

using namespace math;

/**
 * provide comparision operators for the blitz TinyVector.
 * they should be implemented by blitz, but i cannot get
 * them to work?!
 */
bool operator==(math::Vec &t1, math::Vec &t2)
{
  bool b = true;
  for(int i=0; i<3; ++i)
    if (t1(i) != t2(i)) b = false;
  return b;
}

/**
 * != operator
 */
bool operator!=(math::Vec &t1, math::Vec &t2)
{
  return !(t1 == t2);
}

/**
 * Testing leap-frog
 */
int leap_frog_test()
{
  int result = 0;
  
  const int SIZE = 10;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.resize(SIZE);

  // have to add real atoms...
  std::set<int> ex;
  for(int i=0; i<SIZE; ++i)
    the_topology.add_solute_atom("HI", 0, 0, 1.0, 0.0, false, ex, ex);
  
  // initialize everything with zero
  math::Vec t(0.0, 0.0, 0.0), t2(1, 1, 1);

  the_system.pos() = 0.0;
  the_system.vel() = 0.0;
  the_system.force() = 0.0;
  
  // the_topology.mass() = 1.0;

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  simulation_type the_simulation(the_topology, the_system);

  // i need an empty forcefield
  interaction::forcefield<simulation_type> the_forcefield;

  // do a step
  algorithm::leap_frog<simulation_type>
    ::step(the_simulation, the_forcefield, 1.0);

  // check that
  bool b=true;

  for(int i=0; i<SIZE; ++i)
    if (!(the_system.pos()(i) == t)) b = false;

  if (!b) ++result;

  // set initial positions
  for(int i=0; i<SIZE; ++i)
    for(int j=0; j<3; ++j)
      the_system.pos()(i)(j) = i*3+j;
  
  // do a step
  algorithm::leap_frog<simulation_type>
    ::step(the_simulation, the_forcefield, 1.0);

  // check
  for(int i=0; i<SIZE; ++i)
    for(int j=0; j<3; ++j)
      if(the_system.pos()(i)(j) != i*3+j) ++result;

  // add some velocities
  for(int i=0; i<SIZE; ++i)
    for(int j=0; j<3; ++j)
      the_system.vel()(i)(j) = i;

  // do a step
  algorithm::leap_frog<simulation_type>
    ::step(the_simulation, the_forcefield, 1.0);

  // check
  for(int i=0; i<SIZE; ++i)
    for(int j=0; j<3; ++j)
      if (the_system.pos()(i)(j) != i*3+j + i*1.0) ++result;
  
  // set a constant force
  // this only works as long as
  // no interactions are calculated
  // (and the force is not reset to 0.0)
  the_system.force() = 1.0;
  the_system.exchange_force();
  the_system.force() = 1.0;

  the_system.pos() = 0.0;
  the_system.vel() = 0.0;

  // do some steps
  for(int i=0; i<5; ++i)
    algorithm::leap_frog<simulation_type>
      ::step(the_simulation, the_forcefield, 1.0);

  // and back
  // this means i must get rid of the last positions!
  // (because the velocities do not match...)
  the_system.exchange_pos();
  for(int i=0; i<4; ++i)
    algorithm::leap_frog<simulation_type>
      ::step(the_simulation, the_forcefield, -1.0);
  
  // check
  for(int i=0; i<SIZE; ++i)
    for(int j=0; j<3; ++j)
      if(the_system.pos()(i)(j) != 0.0) ++result;


  return result;
}

/**
 * Testing Runge-Kutta
 */
int runge_kutta_test()
{
  int result = 0;
  
  const int SIZE = 5;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.resize(SIZE);

  // have to add real atoms...
  std::set<int> ex;
  for(int i=0; i<SIZE; ++i)
    the_topology.add_solute_atom("HI", 0, 0, 1.0, 0.0, false, ex, ex);
  

  // constant velocity
  the_system.pos() = tensor::i;
  the_system.vel() = tensor::i;
  the_system.force() = 0.0;
  
  // the_topology.mass() = 1.0;

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  simulation_type the_simulation(the_topology, the_system);

  // i need an empty forcefield
  interaction::forcefield<simulation_type> the_forcefield;

  // do a few steps step
  for(int i=0; i<10; ++i)
    algorithm::leap_frog<simulation_type>
      ::step(the_simulation, the_forcefield, 1.0);
  
  // save final position and velocities
  VArray p_fin(the_system.pos());
  VArray v_fin(the_system.vel());

  // reinitialize
  the_system.pos() = tensor::i;
  the_system.vel() = tensor::i;
  the_system.force() = 0.0;
  
  algorithm::runge_kutta<simulation_type> the_step;
  
  for(int i=0; i<10; ++i)
    the_step.step(the_simulation, the_forcefield, 1.0);
  
  // check whether leap-frog and runge-kutta
  // give the same result for constant velocity
  // system.

  bool b = true;
  for(int i=0; i<SIZE; ++i)
    if (p_fin(i) != the_system.pos()(i)) b = false;

  if (!b) ++result;

  return result;
}

int main()
{
  int r1;
  if ((r1 = leap_frog_test())){
    std::cout << "leap_frog_test failed" << std::endl;
  }
  int r2;
  if ((r2 = runge_kutta_test())){
    std::cout << "runge_kutta_test failed" << std::endl;
  }
  
  return r1 + r2;
}

