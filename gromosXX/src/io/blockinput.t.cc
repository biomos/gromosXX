/**
 * @file blockinput.t.cc
 * test routines for blockinput.
 */

#include <iostream>
#include <iomanip>

#include "../math/gmath.h"
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
