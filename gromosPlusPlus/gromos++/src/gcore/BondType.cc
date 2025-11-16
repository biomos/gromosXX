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

// gcore_BondType.cc
#include "BondType.h"

#include <new>

using gcore::BondType;

BondType &BondType::operator=(const BondType &b) {
  if (this != &b) {
    this->~BondType();
    new(this) BondType(b);
  }
  return *this;
}

BondType::BondType(int c, double fc, double l, bool quartic) : d_code(c), d_b0(l) {
  if (quartic) {
    d_fc=fc;
    d_hfc = 2.0 * d_b0 * d_b0 * d_fc;
  } else  {
    d_hfc=fc;
    d_fc = d_hfc / (2 * d_b0 * d_b0);
  }
}

