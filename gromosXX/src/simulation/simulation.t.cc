/**
 * @file simulation.t.cc
 * test routines for the simulation library
 */

#include <iostream>
#include <iomanip>

#include "../math/gmath.h"

#include "simulation/simulation.h"
#include "topology/topology.h"
#include "system/system.h"

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>


using namespace boost::unit_test_framework;

/**
 * some simple tests on the
 * simulation, system and topology classes
 */
void construct_test()
{
  const int SIZE = 10;
  const int SIZE2 = 20;
  const int SIZE3 = 5;
  
  // this is rather a compile test
  simulation::system the_system;
  simulation::topology the_topology;
  
  simulation::simulation<simulation::topology, simulation::system>
    the_simulation(the_topology, the_system);
  
  // two pseudo tests (-> but they failed the first time!!!)
  BOOST_CHECK_EQUAL(&the_system, &the_simulation.system());
  BOOST_CHECK_EQUAL(&the_topology, &the_simulation.topology());

  // resize the system (initially it's 0 size)
  the_simulation.system().resize(SIZE);
  for(int i=0; i<SIZE; ++i)
    the_simulation.system().pos()(i)(0) = i;
  
  // check that
  for(int i=0; i<SIZE; ++i){
    BOOST_CHECK_EQUAL(the_system.pos()(i)(0), i);
  }
  
  // resize again (enlarge)
  the_simulation.system().resize(SIZE2);
  BOOST_CHECK_EQUAL(the_system.pos()(9)(0), 9);
  
  // resize again (shrink)
  the_simulation.system().resize(SIZE3);
  BOOST_CHECK_EQUAL(the_system.pos()(4)(0), 4);

  // does this fail? yes! is it an exception?
  // the_simulation.system().pos()(SIZE2)(0) = 5;

}

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test= BOOST_TEST_SUITE( "simulation test" );
    test->add( BOOST_TEST_CASE( &construct_test ) );

    return test;
}

