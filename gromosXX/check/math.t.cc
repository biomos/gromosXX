/**
 * @file math.t.cc
 * check the simple things in math.
 */

int math_periodicity()
{
  int result = 0;
  int last_result;
  
  CHECKING("math::periodicity<vacuum>", last_result);
  
  Box box(Vec(3.0, 0.0, 0.0),
	     Vec(1.0, 1.0, 0.0),
	     Vec(1.0, 1.0, 2.0));
  
  Periodicity<vacuum> pv;
  pv.box(box);

  Vec v1(0.0, 0.0, 0.0);
  Vec v2(5.0, 2.0, 2.0);
  
  Vec r1;
  pv.nearest_image(v2, v1, r1);
  
  
  Periodicity<math::any> pav(vacuum);
  pav.box(box);
  
  Vec r2;
  pav.nearest_image(v2, v1, r2);

  CHECK_EQUAL(v2, r1, last_result);
  CHECK_EQUAL(r1, r2, last_result);

  RESULT(last_result, result);

  CHECKING("math::periodicity<triclinic>", last_result);

  Periodicity<triclinic> pt;
  pt.box(box);
  
  Vec r3;
  pt.nearest_image(v2, v1, r3);

  Periodicity<math::any> pat(triclinic);
  pat.box(box);
  
  Vec r4;
  pat.nearest_image(v2, v1, r4);
  
  CHECK_EQUAL(v1, r3, last_result);
  CHECK_EQUAL(r3, r4, last_result);

  RESULT(last_result, result);
  
  CHECKING("math::periodicity<triclinic> general vector", last_result);
  Vec v3(1.5, -0.75, -2.0);
  Vec v4(-1.0, -0.25, -0.8);
  
  Vec r5, r6;
  pt.nearest_image(v3, v4, r5);
  pat.nearest_image(v3, v4, r6);

  Vec r7(0.5, 0.5, 0.8);

  CHECK_EQUAL(r5, r7, last_result);
  CHECK_EQUAL(r5, r6, last_result);

  RESULT(last_result, result);

  return result;
}
