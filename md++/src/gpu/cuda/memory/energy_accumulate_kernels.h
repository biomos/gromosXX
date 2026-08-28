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
 * @file energy_accumulate_kernels.h
 * atomicAdd's a GPU-native algorithm's own private per-energy-group
 * buffer into one of gpu::Configuration's shared energy_* device
 * buffers (gpu::MIRROR_ENERGY) -- same rationale and shape as
 * virial_accumulate_kernels.h, just sized to num_energy_groups (a
 * runtime value) instead of a fixed 9. Replaces each bonded/special
 * term's own "cudaStreamSynchronize() + host merge loop reading a
 * gpu::cuvector (managed memory)" pattern, which triggered a real
 * unified-memory page-fault migration every step for what should be a
 * handful of bytes (see git history for the profiling that found
 * this) -- one atomicAdd kernel per contributor, one cudaMemcpy for
 * the whole run's worth of contributors at Energy_Calculation, no
 * per-algorithm host touch at all.
 *
 * `dst` is one of gpu::Configuration's energy_{bond,angle,improper,
 * dihedral,posrest} buffers; `src` is always a plain device-only
 * `double*` scratch buffer the caller wrote via its own kernel, sized
 * to the same num_energy_groups.
 */

#pragma once

namespace gpu {

  void launch_accumulate_energy(
      double* dst, const double* src, unsigned num_groups, cudaStream_t stream = 0);

} // namespace gpu
