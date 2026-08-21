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
 * @file sparse_force_accumulate.h
 * Shared by every bonded-term GPU `Interaction` (quartic bond, angle,
 * improper dihedral, dihedral): downloads the device force buffer into
 * `conf.current().force` for only the atoms a term list actually
 * references, not the whole system.
 *
 * A bonded term's atom set is a tiny fraction of a typical system's
 * total atom count -- bulk solvent in particular is never referenced
 * by any bonded term -- so looping the accumulate over `[0, num_atoms)`
 * (the device force buffer is zeroed there too, but that's a cheap
 * `cudaMemset`, not the concern) means the host pays for a bounds-
 * checked accessor call and a `math::Vec` construction for every atom
 * in the system, the overwhelming majority of which add exactly zero.
 * `build_touched_atoms()` collects the term list's actual atom indices
 * once in `init()`; `accumulate_sparse_forces()` then only visits those
 * every `calculate_interactions()` call.
 */

#pragma once

#include <algorithm>
#include <set>
#include <vector>

#include "gpu/cuda/memory/precision.h"

namespace gpu {

  /**
   * Sorted, de-duplicated atom indices referenced by `terms`, using
   * `atoms_of(term)` to extract each term's own atom indices (e.g.
   * `[](auto const& t){ return std::array{t.i, t.j}; }` for a distance
   * constraint, or `{t.i, t.j, t.k, t.l}` for a dihedral).
   */
  template <typename TermList, typename AtomsOf>
  std::vector<unsigned> build_touched_atoms(const TermList & terms, AtomsOf atoms_of) {
    std::set<unsigned> touched;
    for (const auto & term : terms) {
      for (unsigned a : atoms_of(term)) touched.insert(a);
    }
    return std::vector<unsigned>(touched.begin(), touched.end());
  }

  /**
   * Accumulates (`+=`, matching Forcefield::calculate_interactions()'s
   * zero-once-then-accumulate convention) `force[i]` into
   * `conf.current().force(i)` for every `i` in `touched_atoms`.
   */
  template <typename ConfigurationT>
  void accumulate_sparse_forces(ConfigurationT & conf,
                                 const FPH3_TYPE* force,
                                 const std::vector<unsigned> & touched_atoms) {
    for (unsigned i : touched_atoms) {
      conf.current().force(i) += math::Vec(force[i].x, force[i].y, force[i].z);
    }
  }

} // namespace gpu
