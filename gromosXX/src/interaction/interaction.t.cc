/**
 * @file interaction.t.cc
 * test routines for the interaction library
 */

#include <iostream>
#include <iomanip>

#include "../math/gmath.h"

#include "../simulation/simulation/simulation.h"
#include "../simulation/topology/topology.h"
#include "../simulation/system/system.h"

#include "interaction/interaction.h"
#include "pairlist/simple_pairlist.h"
#include "interaction/nonbonded_interaction.h"


#include "forcefield/forcefield.h"

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
 * Testing the forcefield
 */
void forcefield_test()
{
  const int SIZE = 10;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.num_solute_atoms(SIZE);

  // initialize everything with zero
  Vec t(0.0, 0.0, 0.0), t2(1, 1, 1);

  the_system.pos() = tensor::i;
  the_system.vel() = 0.0;
  the_system.force() = 0.0;
  
  the_topology.mass() = 1.0;

  simulation::simulation<simulation::topology, simulation::system>
    the_simulation(the_topology, the_system);

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  interaction::simple_pairlist<simulation_type> the_pairlist;
  the_pairlist.make_pairlist(the_simulation);
  
  // the_pairlist.print_pairlist(cout);
  
  interaction::simple_pairlist<simulation_type>::iterator it =
    the_pairlist.begin();
  
  /*
  for( ; !it.eol(); ++it){
    std::cout << std::setw(6) << it.i() << setw(6) << it.j() << std::endl;
  }
  */
  
  // let's create a forcefield...

  interaction::forcefield<simulation_type> the_forcefield;

  interaction::nonbonded_interaction<simulation_type> *the_nb_interaction =
    new interaction::nonbonded_interaction<simulation_type>;
  
  the_forcefield.add_interaction(the_nb_interaction);

  // and calculate the forces...
  the_forcefield.calculate_interactions(the_simulation);

  // total force should be 0
  Vec v = sum(the_simulation.system().force());

  bool b = true;
  if (v == t) b=true;
  else b = false;

  BOOST_CHECK_EQUAL(b, true);

}

test_suite*
init_unit_test_suite(int argc, char*argv[])
{
  test_suite* test=BOOST_TEST_SUITE("interaction test");
  test->add(BOOST_TEST_CASE( &forcefield_test));
  
  return test;
}
