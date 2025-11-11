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

/**
 * @file temperature.cc
 * 
 * Implementation of Temperature
 */
#include "Temperature.h"

#include <cassert>
#include <cstdio>

#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gmath/Vec.h"
#include "../gmath/Physics.h"

utils::Temperature::Temperature (const AtomSpecifier &as, double dof)
    : m_as(as), dof(dof)
{
    //std::cout << "# Degree of freedom: = " << dof 
    //          << " Boltzmann = " << gmath::physConst.get_boltzmann() <<std::endl;
}

double
utils::Temperature::temperature(const gcore::System &sys){
    double e_kin = 0.0;
    int num_atoms = m_as.size();
    for (int i = 0; i < num_atoms; i++){
        double mass = m_as.mass(i);
        gmath::Vec vel = m_as.vel(i);
        e_kin += mass * vel.abs2();
    }
    e_kin /= (dof * gmath::physConst.get_boltzmann());
    //e_kin *= 0.5;
    return e_kin;
}
