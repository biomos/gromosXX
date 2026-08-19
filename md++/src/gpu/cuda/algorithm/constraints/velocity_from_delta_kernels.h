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
 * @file velocity_from_delta_kernels.h
 * Shared kernel for `vel(i) = (pos(i) - old_pos(i)) / dt`, the
 * constrained-velocity-correction step every constraint algorithm
 * (SHAKE/M-SHAKE/LINCS/SETTLE) does identically after correcting
 * positions (matches each CPU reference's own tail exactly). Factored
 * out so constraint algorithms migrated onto the GPU-resident mirror
 * (sim.cuda().configuration_view()/mark_gpu_dirty()) can compute this
 * on-device too, instead of downloading positions to host just to
 * recompute a value that's cheaper to do on the GPU where the
 * corrected positions already live.
 *
 * Takes an explicit atom-index list (uploaded once in each caller's
 * init(), from that algorithm's own constrained_atoms()) rather than a
 * dense [first, last) range, since LINCS's constrained atoms are a
 * possibly-non-contiguous set (a constrained chain), unlike M-SHAKE's
 * contiguous solvent range -- one kernel covers both shapes.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"

namespace gpu {

  void launch_velocity_from_delta(
      const FPL3_TYPE* pos,
      const FPL3_TYPE* old_pos,
      FPL3_TYPE* vel,
      const unsigned* atom_indices,
      unsigned num_indices,
      FPL_TYPE dt_i,
      cudaStream_t stream = 0);

} // namespace gpu
