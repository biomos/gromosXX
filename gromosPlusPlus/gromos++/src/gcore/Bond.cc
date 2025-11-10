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

// gcore_Bond.cc

#include "Bond.h"

#include <iostream>

using gcore::Bond;

Bond::Bond(int a, int b, bool warn){
  if(a<b){
    d_a[0]=a;
    d_a[1]=b;
  }
  else{
    d_a[0]=b;
    d_a[1]=a;
    if(warn) {
      std::cerr << "NOTE: order of atoms changed in bond:\n";
      std::cerr << "      " << a+1 << "," << b+1 << " -> " << b+1 << "," << a+1 << std::endl;
    }
  }
  d_type=-1;
}

Bond::Bond(const Bond &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_type=a.d_type;
}

Bond &Bond::operator=(const Bond &b){
  if(this != &b){
    this->~Bond();
    new(this) Bond(b);
  }
  return *this;
}

int gcore::operator<(const Bond &a, const Bond &b){
  if (a[0]<b[0])return 1;
  else if((a[0]==b[0])&&(a[1]<b[1]))return 1;
  else if((a[0]==b[0])&&(a[1]==b[1])&&(a.type()<b.type()))return 1;
  return 0;
}
