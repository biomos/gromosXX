/**
 * @file math.t.cc
 * test routines for the math library
 * (and an example of how to use blitz++)
 */

#include <iostream>
#include <iomanip>

#include "gmath.h"

using namespace math;

/**
 * testing math thingies
 */
int math_test()
{
  int result = 0;
  
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

  return result;
}

int main()
{
  int r1;
  if ((r1 = math_test()))
    std::cout << "math_test failed" << std::endl;
  
  return r1;

}

