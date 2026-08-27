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
 * @file molecular_virial_kernels.h
 * Host-callable entry points for the two-pass GPU-native molecular-
 * virial correction (recovering a molecular virial from the atomic one
 * that bonded/nonbonded/constraint GPU kernels accumulate). See
 * util/prepare_virial.cc's _atomic_to_molecular_virial() for the exact
 * CPU reference this reproduces.
 *
 * Pass 1 (launch_group_com): one thread per pressure group
 * (topo.pressure_groups(), a CSR offsets array -- group g owns atoms
 * [group_offsets[g], group_offsets[g+1]), contiguous by construction).
 * Each thread does the same serial "unwrap and mass-average" walk the
 * CPU's _centre_of_mass() does: prev starts at the first atom's own
 * position; for each atom a in the group, p = nearest_image(pos(a),
 * prev), com_pos accumulates mass(a)*(p+prev), prev += p. This is a
 * genuine sequential dependency within a group (each atom's unwrapped
 * position depends on the previous atom's), not something to further
 * parallelize -- pressure groups are small (individual molecules), so
 * one thread per group is the natural granularity, exactly mirroring
 * how the CPU reference is itself an inherently serial per-group walk.
 *
 * Pass 2 (launch_molecular_virial_correction): one thread per atom.
 * r = nearest_image(pos(a), com_pos[group_id[a]]),
 * corrP(b, a_component) -= force(a)(a_component) * r(b_component)
 * (note the CPU code's `corrP(b,a) += force(a)(a)*r(b)` immediately
 * followed by `virial_tensor -= corrP` -- this kernel accumulates the
 * already-negated contribution directly, so the caller can feed the
 * result straight into launch_accumulate_virial9(), which only adds).
 * group_id[a] is precomputed once at init() (static for a normal run,
 * same as every other per-atom lookup table the bonded-term classes
 * upload once), avoiding a per-kernel-call binary search over group
 * boundaries.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_group_com(
      math::CuVArray::View pos,
      const float* mass,
      const unsigned* group_offsets,
      unsigned num_groups,
      math::boundary_enum boundary,
      math::Box box,
      FPL3_TYPE* com_pos,
      cudaStream_t stream = 0);

  void launch_molecular_virial_correction(
      math::CuVArray::View pos,
      const FPH3_TYPE* force,
      const unsigned* group_id,
      const FPL3_TYPE* com_pos,
      unsigned num_atoms,
      math::boundary_enum boundary,
      math::Box box,
      double* virial,
      cudaStream_t stream = 0);

} // namespace gpu
