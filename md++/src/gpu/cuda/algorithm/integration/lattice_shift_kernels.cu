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
 * @file lattice_shift_kernels.cu
 * GPU-native lattice-shift tracking kernel -- one thread per
 * chargegroup (not per atom): each thread computes its chargegroup's
 * representative position (centre of geometry for solute, first atom
 * for solvent -- exactly math::Periodicity<b>::
 * put_chargegroups_into_box_saving_shifts()'s split, math/periodicity.cc),
 * box-wraps it, then applies + records that same translation for every
 * atom in the group via a small inner loop (solvent chargegroups are
 * typically 1-4 atoms; solute chargegroups similar) -- same shape as
 * gpu::Periodicity<B>::prepare_chargegroup() (gpu/cuda/math/
 * periodicity.h), which the pairlist already uses for its own (shift-
 * unaware) box-wrap.
 *
 * Only gpu::Periodicity<B>::put_into_box() is reused directly; the
 * cog/shift-recording logic is self-contained here (not routed through
 * prepare_chargegroup(), which has no shift-array parameter and is
 * pairlist-owned -- see PLAN.md's "don't extend working perf-critical
 * code for an unrelated feature" convention).
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "lattice_shift_kernels.h"

namespace {

template <math::boundary_enum B>
__global__ void lattice_shift_kernel(
    math::CuVArray::View pos,
    math::CuVArray::View shift,
    const int* __restrict__ cg_offsets,
    unsigned num_chargegroups,
    unsigned num_solute_chargegroups,
    FPL9_TYPE cartesian_to_oblique,
    gpu::Periodicity<B> periodicity)
{
    const unsigned cg_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (cg_i >= num_chargegroups) return;

    const unsigned cg_begin = static_cast<unsigned>(cg_offsets[cg_i]);
    const unsigned cg_end   = static_cast<unsigned>(cg_offsets[cg_i + 1]);

    FPL3_TYPE cog;
    if (cg_i < num_solute_chargegroups) {
        cog = FPL3_TYPE{0, 0, 0};
        for (unsigned a = cg_begin; a < cg_end; ++a) {
            const FPL3_TYPE p = pos(a);
            cog.x += p.x; cog.y += p.y; cog.z += p.z;
        }
        const FPL_TYPE n = static_cast<FPL_TYPE>(cg_end - cg_begin);
        cog.x /= n; cog.y /= n; cog.z /= n;
    } else {
        cog = pos(cg_begin);
    }

    FPL3_TYPE v_box = cog;
    periodicity.put_into_box(v_box);

    const FPL3_TYPE trans = FPL3_TYPE{v_box.x - cog.x, v_box.y - cog.y, v_box.z - cog.z};
    const FPL3_TYPE trans_shift = FPL3_TYPE{
        cartesian_to_oblique(0, 0) * trans.x + cartesian_to_oblique(0, 1) * trans.y + cartesian_to_oblique(0, 2) * trans.z,
        cartesian_to_oblique(1, 0) * trans.x + cartesian_to_oblique(1, 1) * trans.y + cartesian_to_oblique(1, 2) * trans.z,
        cartesian_to_oblique(2, 0) * trans.x + cartesian_to_oblique(2, 1) * trans.y + cartesian_to_oblique(2, 2) * trans.z
    };

    for (unsigned a = cg_begin; a < cg_end; ++a) {
        FPL3_TYPE p = pos(a);
        p.x += trans.x; p.y += trans.y; p.z += trans.z;
        pos(a) = p;

        FPL3_TYPE s = shift(a);
        s.x += trans_shift.x; s.y += trans_shift.y; s.z += trans_shift.z;
        shift(a) = s;
    }
}

} // namespace

template <math::boundary_enum B>
void gpu::launch_lattice_shift(math::CuVArray::View pos,
                                math::CuVArray::View shift,
                                const int* cg_offsets,
                                unsigned num_chargegroups,
                                unsigned num_solute_chargegroups,
                                const FPL9_TYPE & cartesian_to_oblique,
                                const math::Box & box,
                                cudaStream_t stream) {
    const unsigned threads = 256;
    const unsigned blocks  = (num_chargegroups + threads - 1) / threads;
    gpu::Periodicity<B> periodicity(box);
    lattice_shift_kernel<B><<<blocks, threads, 0, stream>>>(
        pos, shift, cg_offsets, num_chargegroups, num_solute_chargegroups,
        cartesian_to_oblique, periodicity);
}

template void gpu::launch_lattice_shift<math::vacuum>(
    math::CuVArray::View, math::CuVArray::View, const int*, unsigned, unsigned,
    const FPL9_TYPE &, const math::Box &, cudaStream_t);
template void gpu::launch_lattice_shift<math::rectangular>(
    math::CuVArray::View, math::CuVArray::View, const int*, unsigned, unsigned,
    const FPL9_TYPE &, const math::Box &, cudaStream_t);
template void gpu::launch_lattice_shift<math::triclinic>(
    math::CuVArray::View, math::CuVArray::View, const int*, unsigned, unsigned,
    const FPL9_TYPE &, const math::Box &, cudaStream_t);
