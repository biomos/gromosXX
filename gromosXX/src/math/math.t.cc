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

int periodicity_test()
{
  int result = 0;
  
  Matrix box(Vec(10.0, 0.0, 0.0),
	     Vec(0.0, 5.0, 0.0),
	     Vec(0.0, 0.0, 5.0));
  
  Periodicity<vacuum> pv(box);
  
  Vec v1(1.0, 2.5, 2.5);
  Vec v2(9.0, 2.5, 2.5);
  Vec v3;
  
  pv.nearest_image(v1, v2, v3);
  
  // std::cout << "nearest image(vacuum): " << v3 << std::endl;

  Vec rv(-8.0, 0.0, 0.0);
  if (v3 != rv) ++result;
  
  Periodicity<triclinic> pi(box);
  pi.nearest_image(v1, v2, v3);

  // std::cout << "nearest image(triclinic): " << v3 << std::endl;
  
  Vec ri(2.0, 0.0, 0.0);
  if (v3 != ri) ++result;

  return 0;
}


int main()
{
  int r1;
  if ((r1 = math_test()))
    std::cout << "math_test failed" << std::endl;
  int r2;
  if ((r2 = periodicity_test()))
    std::cout << "periodicity_test failed" << std::endl;
  
  return r1 + r2;

}

