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
 * @file virial_accumulate_kernels.h
 * atomicAdd's a GPU-native algorithm's own private 9-element virial
 * buffer into the shared GPU mirror's virial_tensor field
 * (gpu::MIRROR_VIRIAL), entirely on-device -- replaces the per-
 * algorithm CPU merge loop (`conf.current()/old().virial_tensor(b,a)
 * += m_virial[...]`) every GPU-native force/constraint algorithm
 * previously used, which required a host round trip every step.
 *
 * Unlike constraint_force_publish_kernels.h (a plain write -- disjoint
 * atom ranges per algorithm), virial_tensor is a single global
 * accumulator every bonded term, NonBonded, *and* the active
 * constraint algorithms all contribute to in the same step, so this
 * genuinely needs atomicAdd, not a plain write. Safe regardless of
 * which stream each caller runs on: the mirror's virial_tensor is
 * zeroed once per step by CudaManager::zero_mirror_force(), before any
 * contributor runs.
 *
 * `dst` is a raw FPH_TYPE* over the mirror's FPH9_TYPE (FP9<FPH_TYPE>,
 * alignas(16), a union of named fields and a 3x3 m[][] array -- see
 * gpu/cuda/memory/float9.h), reinterpreted as a flat 9-element array;
 * `src` is always a plain `double`-typed private buffer (every GPU-
 * native algorithm's own m_virial, gpu::cuvector<double>, sized 9,
 * row-major b*3+a -- SHAKE/SETTLE/M-SHAKE/every bonded term/NonBonded
 * all use this exact shape).
 */

#pragma once

#include "gpu/cuda/memory/precision.h"

namespace gpu {

  void launch_accumulate_virial9(
      FPH_TYPE* dst, const double* src, cudaStream_t stream = 0);

} // namespace gpu
