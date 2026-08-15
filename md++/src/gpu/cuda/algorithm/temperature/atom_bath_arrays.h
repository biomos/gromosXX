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
 * @file atom_bath_arrays.h
 * Per-atom temperature-group index and (com_bath, ir_bath) indices,
 * shared by every `Thermostat`-derived `gpuBackend` specialization that
 * needs `Thermostat::scale()`'s uniform per-atom formula
 * (berendsen_thermostat_gpu.cc, nosehoover_thermostat_gpu.cc) --
 * both feed the exact same `launch_thermostat_scale_apply` kernel, so
 * they need the exact same per-atom lookup tables.
 *
 * Built once from `topo.temperature_groups()` and
 * `sim.multibath().in_bath()` (the existing CPU-side lookup, reused
 * rather than re-deriving its range logic). Cached in a function-local
 * static, same rationale as Remove_COM_Motion's/Temperature_
 * Calculation's scratch buffers: topology and multibath structure are
 * static for a normal run. Since at most one thermostat is ever active
 * in a given run (create_md_sequence.cc picks exactly one), the two
 * call sites sharing this single cache (an `inline` function's local
 * statics have one definition across the whole program) costs nothing
 * -- it simply avoids rebuilding the same tables twice if a future
 * caller ever needed both.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"

namespace topology { class Topology; }
namespace simulation { class Simulation; }

namespace gpu {

  struct AtomBathArrays {
    gpu::cuvector<unsigned> group_index;
    gpu::cuvector<unsigned> com_bath;
    gpu::cuvector<unsigned> ir_bath;
    unsigned num_groups = 0;
  };

  inline const AtomBathArrays & atom_bath_arrays(const topology::Topology & topo,
                                                  const simulation::Simulation & sim) {
    static AtomBathArrays arrays;
    static unsigned cached_num_atoms = 0;
    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    if (cached_num_atoms == num_atoms) return arrays;

    arrays.group_index.resize(num_atoms);
    arrays.com_bath.resize(num_atoms);
    arrays.ir_bath.resize(num_atoms);

    // temperature_groups() has num_groups+1 entries with a leading 0
    // and EXCLUSIVE end boundaries (see temperature_calculation_gpu.cc's
    // atom_temp_group() for the confirmation from in_topology.cc) --
    // NOT the same convention as energy_groups(). Group g spans atoms
    // [tg[g], tg[g+1]).
    const std::vector<unsigned int> & tg = topo.temperature_groups();
    unsigned group = 0;
    for (unsigned atom = 0; atom < num_atoms; ++atom) {
      if (group + 1 < tg.size() && atom == tg[group + 1]) ++group;
      arrays.group_index[atom] = group;
      unsigned com = 0, ir = 0;
      sim.multibath().in_bath(atom, com, ir);
      arrays.com_bath[atom] = com;
      arrays.ir_bath[atom]  = ir;
    }
    arrays.num_groups = static_cast<unsigned>(tg.size()) - 1;
    cached_num_atoms = num_atoms;
    return arrays;
  }

} // namespace gpu
