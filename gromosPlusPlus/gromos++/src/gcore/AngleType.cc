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

// gcore_AngleType.cc
#include "AngleType.h"

#include <iostream>
#include <cmath>
#include <sstream>

#include "GromosForceField.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../gmath/Physics.h"

using namespace std;
using gcore::AngleType;

AngleType &AngleType::operator=(const AngleType &b) {
  if (this != &b) {
    this->~AngleType();
    new(this) AngleType(b);
  }
  return *this;
}

AngleType::AngleType(int c, double fc, double l) : d_code(c), d_t0(l), d_fc(fc),
d_afc(0.0) {

  // first let's make sure it is not simply zero.
  if (d_fc == 0.0) {
    return;
  }

  double t_0 = d_t0 * gmath::physConst.get_degree2radian();

  // this is kT at 300K - a reasonable choice for biomolecular simulations.
  const double kT = gmath::physConst.get_boltzmann() * 300.0;

  // the two solutions for theta where V_harm(theta) = V_cos-harm(theta)
  if(d_fc < 0) {
    stringstream msg;
    msg << "negative cosine-harmonic force constant k = " << d_fc
            << " does not make sense";
    throw gromos::Exception("AngleType", msg.str());
  }
  
  double costheta1 = cos(t_0) + sqrt(kT / d_fc);
  double costheta2 = cos(t_0) - sqrt(kT / d_fc);
  
  double theta1, theta2;
  
  if(costheta1 <= 1 && costheta1 >= -1) {
    theta1 = acos(costheta1);
    if(costheta2 <= 1 && costheta2 >= -1) {
      theta2 = acos(costheta2);
    } else {
      theta2 = 2 * t_0 - theta1;
    }
  } else if (costheta2 <= 1 && costheta2 >= -1) {
    theta2 = acos(costheta2);
    theta1 = 2 * t_0 - theta2;
  } else {
    throw gromos::Exception("AngleType", "Cannot convert cosine-harmonic bond "
            "angle to a pure harmonic one");
  }
  
  double term1 = (theta1 - t_0) * (theta1 - t_0);
  double term2 = (theta2 - t_0) * (theta2 - t_0);
  
  d_afc = 2 * kT / (term1 + term2) * gmath::physConst.get_degree2radian() * gmath::physConst.get_degree2radian();
  
  // do a warning as this may be a bit tricky
  if (!args::Arguments::outG96) {
    std::cerr << "Warning for bond angle type " << d_code + 1 << ": Harmonic force "
            "constant was calculated at 300K." << std::endl;
  }
}
