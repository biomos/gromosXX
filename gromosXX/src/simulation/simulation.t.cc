/**
 * @file simulation.t.cc
 * test routines for the simulation library
 */

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "../math/gmath.h"

#include "simulation.h"

// #include <boost/test/unit_test_suite.hpp>
// #include <boost/test/test_tools.hpp>


// using namespace boost::unit_test_framework;

/**
 * some simple tests on the
 * simulation, system and topology classes
 */
int construct_test()
{
  int result = 0;
  
  const int SIZE = 10;
  const int SIZE2 = 20;
  const int SIZE3 = 5;
  
  // this is rather a compile test
  simulation::system the_system;
  simulation::topology the_topology;
  
  simulation::simulation<simulation::topology, simulation::system>
    the_simulation(the_topology, the_system);
  
  // two pseudo tests (-> but they failed the first time!!!)
  if (&the_system != &the_simulation.system()) ++result;
  if (&the_topology != &the_simulation.topology()) ++result;

  // resize the system (initially it's 0 size)
  the_simulation.system().resize(SIZE);
  for(int i=0; i<SIZE; ++i)
    the_simulation.system().pos()(i)(0) = i;
  
  // check that
  for(int i=0; i<SIZE; ++i){
    if (the_system.pos()(i)(0) != i) ++result;
  }
  
  // resize again (enlarge)
  the_simulation.system().resize(SIZE2);
  if (the_system.pos()(9)(0) != 9) ++result;
  
  // resize again (shrink)
  the_simulation.system().resize(SIZE3);
  if (the_system.pos()(4)(0) != 4) ++result;

  // does this fail? yes! is it an exception?
  // the_simulation.system().pos()(SIZE2)(0) = 5;

  // try to add a bond
  the_simulation.topology().bonds().add(1, 2, 1);
  simulation::bond::iterator b_it = 
    the_simulation.topology().bonds().begin();
  
  for( ; !b_it.eol(); ++b_it)
    if (1 != b_it.i()) ++result;
  
  return result;
}

/*
test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test= BOOST_TEST_SUITE( "simulation test" );
    test->add( BOOST_TEST_CASE( &construct_test ) );

    return test;
}
*/

int main()
{
  int result = 0;
  
  result += construct_test();
  
  return result;
}


