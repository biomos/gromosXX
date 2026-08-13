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

      /**
       * CSR exclusion list: solute atom i (i < num_solute_atoms) excludes
       * atoms excl_list[excl_ptr[i] .. excl_ptr[i+1]), sorted, all > i.
       * Solvent atoms have no entries -- solvent-solvent tiles never need
       * exclusion masking because different solvent chargegroups are
       * different molecules by construction (see TILE_PAIRLIST_DESIGN.md
       * §1, matching how Standard_Pairlist_Algorithm's _solvent_solvent
       * never calls excluded_solute_pair either).
       */
      const int* const excl_ptr;
      const int* const excl_list;
      const unsigned num_solute_atoms;

      HOSTDEVICE TopologyView() = default;

      HOSTDEVICE TopologyView(
          int* iac_, float* mass_, float* imass_, float* charge_, int* chargegroup_,
          unsigned n_s_cg, unsigned n_s_mol, unsigned n_atoms, unsigned n_cg,
          int* excl_ptr_, int* excl_list_, unsigned n_s_atoms)
          : iac(iac_), mass(mass_), inverse_mass(imass_), charge(charge_), chargegroup(chargegroup_),
            num_solute_chargegroups(n_s_cg), num_solute_molecules(n_s_mol),
            num_atoms(n_atoms), num_chargegroups(n_cg),
            excl_ptr(excl_ptr_), excl_list(excl_list_), num_solute_atoms(n_s_atoms) {}

      /**
       * @brief is atom j excluded from atom i's nonbonded interactions?
       * Requires i < j (same convention as CPU excluded_solute_pair).
       */
      HOSTDEVICE bool is_excluded(unsigned i, unsigned j) const {
          if (i >= num_solute_atoms) return false;
          int lo = excl_ptr[i];
          int hi = excl_ptr[i + 1] - 1;
          const int jj = static_cast<int>(j);
          while (lo <= hi) {
              const int mid = lo + (hi - lo) / 2;
              const int v = excl_list[mid];
              if (v == jj) return true;
              if (v < jj) lo = mid + 1;
              else hi = mid - 1;
          }
          return false;
      }
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

    /**
     * CSR exclusion list, see TopologyView. Sized num_solute_atoms+1 /
     * num_exclusion_entries respectively. Built once from
     * topo.all_exclusion(i), i < num_solute_atoms, in the constructor only
     * -- exclusions are static for a normal run, so update() (called for
     * per-step atom-property resync) intentionally does not touch these.
     */
    int*       excl_ptr;
    int*       excl_list;
    unsigned   num_solute_atoms;
    unsigned   num_exclusion_entries;

    Topology(const topology::Topology& topo);
    ~Topology();

    void update(const topology::Topology & topo);

    // Construct a trivial view to pass to kernels
    const View view() const {
        return View{
            iac, mass, inverse_mass, charge, chargegroup,
            num_solute_chargegroups, num_solute_molecules,
            num_atoms, num_chargegroups,
            excl_ptr, excl_list, num_solute_atoms
        };
    }
};


} // namespace configuration