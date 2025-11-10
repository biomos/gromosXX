/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// gmath_Vec.t.cc

#include "Vec.h"

#include <iostream>

using gmath::Vec;
using namespace std;

ostream &operator<<(ostream &o, const Vec &v){
  o << '(' << v[0] << ' '
    << v[1] << ' '
    << v[2] << ')';
  return o;
}

int main(){
  Vec v(1,2,3);
  Vec w(4,5,6);
  cout << "v = "<< v << " w = " << w << endl;
  cout << "v + w = "<< v+w << endl;
  cout << "v - w = "<< v-w << endl;
  cout << "v * 3 = "<< v*3 << endl;
  cout << "3 * v = "<< 3*v << endl;
  cout << "v / 2 = " << v/2 << endl;
  w-=v;
  cout << "w-=v: " << w << endl;
  w+=v;
  cout << "w+=v: " << w << endl;
  cout << "v.dot(w) = " << v.dot(w) << endl;
  cout << "v.cross(w) = " << v.cross(w) << endl;
  Vec t;
  t=w+v;
  cout << "t=w+v: " << t << endl;
  cout << "v.abs2() = " << v.abs2() << endl;
  cout << "v.abs() = " << v.abs() << endl;
  return 0;
}

