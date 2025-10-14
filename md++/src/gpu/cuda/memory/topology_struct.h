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
#include "gpu/cuda/cuhostdevice.h"

namespace topology {
  class Topology;
}

namespace gpu {
  struct TopologyView {
      const int* const iac;
      const float* const mass;
      const float* const inverse_mass;
      const float* const charge;
      const int* const chargegroup;

      const unsigned num_solute_chargegroups;
      const unsigned num_solute_molecules;
      const unsigned num_atoms;
      const unsigned num_chargegroups;

      HOSTDEVICE TopologyView() = default;

      HOSTDEVICE TopologyView(
          int* iac_, float* mass_, float* imass_, float* charge_, int* chargegroup_,
          unsigned n_s_cg, unsigned n_s_mol, unsigned n_atoms, unsigned n_cg)
          : iac(iac_), mass(mass_), inverse_mass(imass_), charge(charge_), chargegroup(chargegroup_),
            num_solute_chargegroups(n_s_cg), num_solute_molecules(n_s_mol),
            num_atoms(n_atoms), num_chargegroups(n_cg) {}
  };

  /**
   * @brief Holds GPU-side topology data.
   */
  struct Topology {
    using View = TopologyView;
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

    // Construct a trivial view to pass to kernels
    const View view() const {
        return View{
            iac, mass, inverse_mass, charge, chargegroup,
            num_solute_chargegroups, num_solute_molecules,
            num_atoms, num_chargegroups
        };
    }
};


} // namespace configuration