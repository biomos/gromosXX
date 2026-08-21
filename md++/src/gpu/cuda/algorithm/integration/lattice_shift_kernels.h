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
 * @file lattice_shift_kernels.h
 * Host-callable entry points for the GPU-native lattice-shift tracking
 * kernel (put chargegroups into the box, record the cumulative shift
 * per atom). One kernel launcher per math::boundary_enum, matching how
 * the pairlist's own boundary dispatch works (SPLIT_BOUNDARY at the
 * host call site, gpu::Periodicity<BOUNDARY> as the device-side
 * template parameter) -- there is no single runtime-polymorphic device
 * function for this, so the host wrapper is instantiated once per
 * boundary type instead of taking an enum parameter.
 *
 * Mirrors math::Periodicity<b>::put_chargegroups_into_box_saving_shifts()
 * (math/periodicity.cc) exactly -- see that function for the derivation
 * (centre-of-geometry for solute chargegroups, first-atom position for
 * solvent, translate into the box, record the translation in oblique
 * (lattice-vector) coordinates).
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @brief Box-wrap every chargegroup's atoms in place (pos) and
   * accumulate the applied translation, expressed in oblique
   * (lattice-vector) coordinates, into `shift` for every atom in that
   * chargegroup. `cg_offsets` is the standard CSR-style chargegroup
   * boundary array (size num_chargegroups+1, chargegroup i spans atoms
   * [cg_offsets[i], cg_offsets[i+1])) -- the same array already
   * uploaded once as part of gpu::Topology::View::chargegroup, no
   * separate upload needed. `cartesian_to_oblique` is the 3x3 matrix
   * computed on the host once per call from the current box (see
   * Lattice_Shift_Tracker<gpuBackend>::apply()) -- O(9) work, not
   * worth a kernel.
   */
  template <math::boundary_enum B>
  void launch_lattice_shift(math::CuVArray::View pos,
                             math::CuVArray::View shift,
                             const int* cg_offsets,
                             unsigned num_chargegroups,
                             unsigned num_solute_chargegroups,
                             const FPL9_TYPE & cartesian_to_oblique,
                             const math::Box & box,
                             cudaStream_t stream = 0);

}
