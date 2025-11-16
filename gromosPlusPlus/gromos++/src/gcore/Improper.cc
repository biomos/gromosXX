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

// gcore_Improper.cc

#include "Improper.h"

#include<iostream>

using gcore::Improper;

Improper::Improper(int a, int b, int c, int d, bool warn){
  if(b<c){
    d_a[0]=a; 
    d_a[1]=b; 
    d_a[2]=c; 
    d_a[3]=d;
  } 
  else{
    d_a[0]=d;
    d_a[1]=c;
    d_a[2]=b;
    d_a[3] = a;
    if (warn) {
      std::cerr << "NOTE: order of atoms changed in improper dihedral:\n";
      std::cerr << "      " << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << " -> "
              << d + 1 << "," << c + 1 << "," << b + 1 << "," << a + 1 << std::endl;
    }
  }
  d_type = -1;
}

Improper::Improper(const Improper &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_a[2]=a.d_a[2];
  d_a[3]=a.d_a[3];
  d_type=a.d_type;
}

Improper &Improper::operator=(const Improper &b){
  if(this != &b){
    this->Improper::~Improper();
    new(this) Improper(b);
  }
  return *this;
}

int gcore::operator<(const Improper &a, const Improper &b){
return (a[0]<b[0]||(a[0]==b[0]&&a[1]<b[1])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]<b[2])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]<b[3])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]==b[3]&&a.type()<b.type()));
}
