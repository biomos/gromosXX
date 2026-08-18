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
 * @file one_four_kernels.cu
 * 1,4-pair ("LJ exception") LJ+CRF kernel. See one_four_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "one_four_kernels.h"

namespace gpu {

/**
 * One thread per solute atom i (grid-stride loop), looping
 * `one_four_list[one_four_ptr[i] .. one_four_ptr[i+1])` -- always j > i,
 * so every 1,4 pair is visited exactly once. Direct global atomics, not
 * shared-memory bucketing -- same reasoning as
 * rf_excluded_solute_kernel: a solute atom's 1,4-pair count is bounded
 * by local bonded topology, not cutoff-sphere atom count.
 */
template <math::boundary_enum BOUNDARY>
__global__ void one_four_kernel(
    const int* __restrict__ one_four_ptr,
    const int* __restrict__ one_four_list,
    unsigned num_solute_atoms,
    math::CuVArray::View pos,
    const int* __restrict__ iac,
    const FPL_TYPE* __restrict__ charge,
    const unsigned* __restrict__ atom_energy_group,
    gpu::LJParamView lj,
    gpu::NbSimParams nb,
    FPL_TYPE coulomb_scaling,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* e_lj_total,
    double* e_crf_total,
    double* virial_total) {

    const unsigned num_groups = nb.num_energy_groups;

    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
         i < num_solute_atoms; i += blockDim.x * gridDim.x) {

        const unsigned eg_i = atom_energy_group[i];
        const FPL_TYPE qi = charge[i];

        FPL3_TYPE f_accum{FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0)};

        for (int k = one_four_ptr[i]; k < one_four_ptr[i + 1]; ++k) {
            const unsigned j = static_cast<unsigned>(one_four_list[k]);

            const FPL3_TYPE rvec  = periodicity.nearest_image(pos(i), pos(j));
            const FPL_TYPE dist2  = abs2(rvec);
            const FPL_TYPE dist2i = FPL_TYPE(1) / dist2;
            const FPL_TYPE dist6i = dist2i * dist2i * dist2i;

            const FPL2_TYPE ljp = lj.get_scaled(iac[i], iac[j]);
            const FPL_TYPE cs6  = ljp.x;
            const FPL_TYPE cs12 = ljp.y;
            const FPL_TYPE cs12_dist6i = cs12 * dist6i;

            const FPL_TYPE e_lj = (cs12_dist6i - cs6) * dist6i;

            const FPL_TYPE disti = sqrtf(dist2i);
            const FPL_TYPE q_eps = qi * charge[j] * nb.four_pi_eps_i;

            const FPL_TYPE e_crf = q_eps * (disti * coulomb_scaling -
                                             nb.crf_2cut3i * dist2 - nb.crf_cut);

            const FPL_TYPE f_scalar =
                (cs12_dist6i + cs12_dist6i - cs6) * FPL_TYPE(6) * dist6i * dist2i +
                q_eps * (disti * coulomb_scaling * dist2i + FPL_TYPE(2) * nb.crf_2cut3i);
            const FPL3_TYPE fr = f_scalar * rvec;

            f_accum.x += fr.x;
            f_accum.y += fr.y;
            f_accum.z += fr.z;
            atomicAdd(&force[j].x, -fr.x);
            atomicAdd(&force[j].y, -fr.y);
            atomicAdd(&force[j].z, -fr.z);

            const unsigned eg_j = atom_energy_group[j];
            const unsigned bucket = eg_i * num_groups + eg_j;
            atomicAdd(&e_lj_total[bucket],  static_cast<double>(e_lj));
            atomicAdd(&e_crf_total[bucket], static_cast<double>(e_crf));

            atomicAdd(&virial_total[0], static_cast<double>(rvec.x * fr.x));
            atomicAdd(&virial_total[1], static_cast<double>(rvec.x * fr.y));
            atomicAdd(&virial_total[2], static_cast<double>(rvec.x * fr.z));
            atomicAdd(&virial_total[3], static_cast<double>(rvec.y * fr.x));
            atomicAdd(&virial_total[4], static_cast<double>(rvec.y * fr.y));
            atomicAdd(&virial_total[5], static_cast<double>(rvec.y * fr.z));
            atomicAdd(&virial_total[6], static_cast<double>(rvec.z * fr.x));
            atomicAdd(&virial_total[7], static_cast<double>(rvec.z * fr.y));
            atomicAdd(&virial_total[8], static_cast<double>(rvec.z * fr.z));
        }

        atomicAdd(&force[i].x, f_accum.x);
        atomicAdd(&force[i].y, f_accum.y);
        atomicAdd(&force[i].z, f_accum.z);
    }
}

} // namespace gpu

void gpu::launch_one_four(
    const int* one_four_ptr,
    const int* one_four_list,
    unsigned num_solute_atoms,
    math::CuVArray::View pos,
    const int* iac,
    const FPL_TYPE* charge,
    const unsigned* atom_energy_group,
    gpu::LJParamView lj,
    gpu::NbSimParams nb,
    FPL_TYPE coulomb_scaling,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* force,
    double* e_lj_total,
    double* e_crf_total,
    double* virial_total,
    cudaStream_t stream) {

    if (num_solute_atoms == 0) return;

    constexpr unsigned kThreads = 128;
    const unsigned blocks = (num_solute_atoms + kThreads - 1) / kThreads;

    switch (boundary) {
        case math::vacuum: {
            const gpu::Periodicity<math::vacuum> periodicity(box);
            gpu::one_four_kernel<math::vacuum><<<blocks, kThreads, 0, stream>>>(
                one_four_ptr, one_four_list, num_solute_atoms,
                pos, iac, charge, atom_energy_group, lj, nb, coulomb_scaling,
                periodicity, force, e_lj_total, e_crf_total, virial_total);
            break;
        }
        case math::rectangular: {
            const gpu::Periodicity<math::rectangular> periodicity(box);
            gpu::one_four_kernel<math::rectangular><<<blocks, kThreads, 0, stream>>>(
                one_four_ptr, one_four_list, num_solute_atoms,
                pos, iac, charge, atom_energy_group, lj, nb, coulomb_scaling,
                periodicity, force, e_lj_total, e_crf_total, virial_total);
            break;
        }
        case math::triclinic: {
            const gpu::Periodicity<math::triclinic> periodicity(box);
            gpu::one_four_kernel<math::triclinic><<<blocks, kThreads, 0, stream>>>(
                one_four_ptr, one_four_list, num_solute_atoms,
                pos, iac, charge, atom_energy_group, lj, nb, coulomb_scaling,
                periodicity, force, e_lj_total, e_crf_total, virial_total);
            break;
        }
        default:
            // Unsupported boundary -- matches CUDA_Pairlist_Algorithm::init()'s
            // vacuum/rectangular-only v1 scope; no kernel launch.
            break;
    }
    CUDA_CHECK_ERROR("one_four_kernel");
}
