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
 * @file topology_struct.h
 * A light-weight struct for GPU holding essential Topology data
 * The struct holds all device pointers
 */

#pragma once

// #include "types.h"
// #include "cuvector.h"

namespace topology {
  class Topology;
}

namespace gpu {
  /**
   * @brief Holds GPU-side view of topology data.
   * 
   * All arrays are in single precision (float or float3) to match GPU performance requirements.
   * We assume constants, so we copy only once in the beginning
   */
  struct Topology {
    void*      memory_block; // base pointer for deallocation
    // Atom codes, masses, charges
    int*       iac;
    float*     mass;
    float*     inverse_mass;
    float*     charge;
    int*       chargegroup;
    unsigned   num_solute_chargegroups;
    unsigned   num_solute_molecules;
    unsigned   num_atoms;
    unsigned   num_chargegroups;

    Topology(const topology::Topology& topo);
    ~Topology();

    void update(const topology::Topology & topo);
  };

} // namespace configuration