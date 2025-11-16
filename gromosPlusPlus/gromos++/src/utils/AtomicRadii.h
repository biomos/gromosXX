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
 * @file AtomicRadii.h
 * AtomicRadii methods
 */

#ifndef INCLUDED_UTILS_ATOMICRADII
#define INCLUDED_UTILS_ATOMICRADII

namespace gcore {
  class System;
  class GromosForceField;
}

namespace utils
{
  /**
   * Compute the atomic radii as the minimal Lennard Jones energy distance of the
   * atoms and the probe particle minus the radius of that particle.
   *
   * @param probe_iac the IAC of the probe
   * @param probe_radius the radius of the probe
   * @param sys the system
   * @param gff the GROMOS force field parameters
   */
  void compute_atomic_radii_vdw(gcore::System & sys, const gcore::GromosForceField & gff);
  void compute_atomic_radii_vdw(int probe_iac, double probe_radius, gcore::System & sys, const gcore::GromosForceField & gff);
}

#endif

