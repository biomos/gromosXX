/**
 * @file math.t.cc
 * test routines for the math library
 * (and an example of how to use blitz++)
 */

#include <iostream>
#include <iomanip>

#include "gmath.h"

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

using namespace boost::unit_test_framework;
using namespace math;

/**
 * testing math thingies
 */
void math_test()
{
  const int SIZE = 10;
  
  // Vector expressions
  Vec v1(1.0 , 2.0, 3.0);
  Vec v2(1.0 , 1.0, 1.0);
  double d;
  d = dot(v1, v2);
  Vec v3 = cross(v1, v2);
  v3 = v1 + v2;
  
  // scalar arrays
  SArray m(5);
  m = 2.5;
  m *= 0.25;
  SArray mm(5);
  mm = 1, 2, 3, 4, 5;
  m += mm;
  
  // vector arrays
  VArray x(5), y(5);
  x = 1; // constant
  y = v1; // 5 v1 vectors
  VArray z(x + y);
}

test_suite*
init_unit_test_suite(int argc, char*argv[])
{
  test_suite* test=BOOST_TEST_SUITE("math test");
  test->add(BOOST_TEST_CASE( &math_test));
  
  return test;
}
