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

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

using namespace boost::unit_test_framework;
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
void leap_frog_test()
{
  const int SIZE = 10;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.resize(SIZE);
  the_topology.num_solute_atoms(SIZE);
  
  // initialize everything with zero
  blitz::TinyVector<double, 3> t(0.0, 0.0, 0.0), t2(1, 1, 1);

  the_system.pos() = 0.0;
  the_system.vel() = 0.0;
  the_system.force() = 0.0;
  
  the_topology.mass() = 1.0;

  simulation::simulation<simulation::topology, simulation::system>
    the_simulation(the_topology, the_system);

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  // i need an empty forcefield
  interaction::forcefield<simulation_type> the_forcefield;

  // do a step
  algorithm::leap_frog<simulation_type>
    ::step(the_simulation, the_forcefield, 1.0);

  // check that
  bool b=true;

  for(int i=0; i<SIZE; ++i)
    if (!(the_system.pos()(i) == t)) b = false;

  //firstIndex i;
  //_bz_bool bb = all((the_system.pos() == t), i);
  

  BOOST_CHECK(b);

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
      BOOST_CHECK_EQUAL(the_system.pos()(i)(j),i*3+j);

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
      BOOST_CHECK_EQUAL(the_system.pos()(i)(j),i*3+j + i*1.0);
  
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
      BOOST_CHECK_EQUAL(the_system.pos()(i)(j), 0.0);

}

/**
 * Testing Runge-Kutta
 */
void runge_kutta_test()
{
  const int SIZE = 5;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.resize(SIZE);
  the_topology.num_solute_atoms(SIZE);

  // constant velocity
  the_system.pos() = tensor::i;
  the_system.vel() = tensor::i;
  the_system.force() = 0.0;
  
  the_topology.mass() = 1.0;

  simulation::simulation<simulation::topology, simulation::system>
    the_simulation(the_topology, the_system);

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

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
  
  // std::cout << p_fin << std::endl;
  // std::cout << the_system.pos() << std::endl;

  // check whether leap-frog and runge-kutta
  // give the same result for constant velocity
  // system.

  bool b = true;
  for(int i=0; i<SIZE; ++i)
    if (p_fin(i) != the_system.pos()(i)) b = false;

  BOOST_CHECK_EQUAL(b, true);

}


test_suite*
init_unit_test_suite(int argc, char*argv[])
{
  test_suite* test=BOOST_TEST_SUITE("algorithm test");
  test->add(BOOST_TEST_CASE( &leap_frog_test));
  test->add(BOOST_TEST_CASE( &runge_kutta_test));
  
  return test;
}
