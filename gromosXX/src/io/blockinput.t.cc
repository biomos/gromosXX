/**
 * @file blockinput.t.cc
 * test routines for blockinput.
 */

#include "../math/gmath.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "../simulation/simulation.h"
#include "../interaction/interaction.h"

#include "blockinput.h"
#include "GInStream.h"
#include "topology/InTopology.h"

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

using namespace boost::unit_test_framework;

/**
 * test the io.
 */
void test_blockio()
{
  simulation::system the_system;
  simulation::topology the_topology;
  
  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  simulation_type the_simulation(the_topology, the_system);

  // i need an empty forcefield
  interaction::forcefield<simulation_type> the_forcefield;

  interaction::harmonic_bond_interaction<simulation_type> *bond_interaction
    = new interaction::harmonic_bond_interaction<simulation_type>;
  
  std::ifstream topo_file("/home/markus/test/hexa10.topo");
  io::InTopology topo(topo_file);
  
  topo >> *bond_interaction;

  topo >> the_topology;

  // no tests...
  BOOST_CHECK_EQUAL(1, 1);
}

test_suite*
init_unit_test_suite(int argc, char* argv[])
{
  test_suite* test=BOOST_TEST_SUITE("blockio test");
  test->add(BOOST_TEST_CASE( &test_blockio));

  return test;
}
