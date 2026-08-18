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
 * @file rf_excluded_kernels.cu
 * Reaction-field-for-excluded-pairs kernels. See rf_excluded_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "rf_excluded_kernels.h"

namespace gpu {

/**
 * @brief One thread per solute atom i (grid-stride loop). Self-term
 * (r=0, energy only, exactly half-counted per atom -- matches
 * `RF_excluded_interaction_innerloop`'s `0.5 * e_crf` into the diagonal
 * [eg_i][eg_i] bucket) plus one real force/energy/virial contribution
 * per entry in `excl_list[excl_ptr[i] .. excl_ptr[i+1])` (always j > i,
 * so every excluded pair is visited exactly once, from i's side).
 *
 * No shared-memory energy/virial bucketing (unlike lj_crf_tile_kernel):
 * a solute atom's exclusion-list length is bounded by its local bonded
 * topology (1-2/1-3/1-4 neighbours), not by cutoff-sphere atom count, so
 * there's no O(n^2)-pairs-per-block case to amortize atomics against --
 * a handful of direct atomicAdds per atom is already cheap.
 */
template <math::boundary_enum BOUNDARY>
__global__ void rf_excluded_solute_kernel(
    const int* __restrict__ excl_ptr,
    const int* __restrict__ excl_list,
    unsigned num_solute_atoms,
    math::CuVArray::View pos,
    const FPL_TYPE* __restrict__ charge,
    const unsigned* __restrict__ atom_energy_group,
    gpu::NbSimParams nb,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* e_crf_total,
    double* virial_total) {

    const unsigned num_groups = nb.num_energy_groups;

    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
         i < num_solute_atoms; i += blockDim.x * gridDim.x) {

        const unsigned eg_i = atom_energy_group[i];
        const FPL_TYPE qi = charge[i];

        // self-term: r = 0, so only the distance-independent part of
        // rf_interaction survives: e_crf = qi*qi * four_pi_eps_i * (-crf_cut).
        const double e_self = 0.5 * static_cast<double>(qi) * static_cast<double>(qi) *
                               static_cast<double>(nb.four_pi_eps_i) *
                               (-static_cast<double>(nb.crf_cut));
        atomicAdd(&e_crf_total[eg_i * num_groups + eg_i], e_self);

        // Accumulate atom i's own force contribution locally across all
        // its exclusions, one atomicAdd at the end instead of one per
        // pair -- atom j's side still needs a per-pair atomicAdd since
        // several different i's exclusion lists can reference the same j.
        FPL3_TYPE f_accum{FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0)};

        for (int k = excl_ptr[i]; k < excl_ptr[i + 1]; ++k) {
            const unsigned j = static_cast<unsigned>(excl_list[k]);

            const FPL3_TYPE rvec = periodicity.nearest_image(pos(i), pos(j));
            const FPL_TYPE dist2 = abs2(rvec);
            const FPL_TYPE q_eps = qi * charge[j] * nb.four_pi_eps_i;

            // rf_interaction (nonbonded_term.cc): e_crf = q_eps*(-crf_2cut3i*r^2 - crf_cut),
            // force = q_eps * crf_cut3i * r, crf_cut3i == 2*crf_2cut3i by construction.
            const FPL_TYPE e_crf = q_eps * (-nb.crf_2cut3i * dist2 - nb.crf_cut);
            const FPL_TYPE f_scalar = q_eps * (FPL_TYPE(2) * nb.crf_2cut3i);
            const FPL3_TYPE fr = f_scalar * rvec;

            f_accum.x += fr.x;
            f_accum.y += fr.y;
            f_accum.z += fr.z;
            atomicAdd(&force[j].x, -fr.x);
            atomicAdd(&force[j].y, -fr.y);
            atomicAdd(&force[j].z, -fr.z);

            const unsigned eg_j = atom_energy_group[j];
            atomicAdd(&e_crf_total[eg_i * num_groups + eg_j], static_cast<double>(e_crf));

            // Atomic virial, same convention/index layout as
            // lj_crf_tile_kernel (virial_total[b*3+a] += r(b)*force(a)).
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

/**
 * @brief One thread per solvent chargegroup (grid-stride loop over
 * [num_solute_chargegroups, num_chargegroups)). A solvent chargegroup is
 * one rigid molecule (few atoms, e.g. 3 for SPC water) -- looping its
 * O(k^2) internal pairs on a single thread is cheap and matches the CPU
 * innerloop's own per-chargegroup granularity exactly; the actual
 * parallelism is across the (typically many thousand) chargegroups, one
 * per thread, not within one.
 *
 * Energy only (rigid solvent, no internal forces -- matches the CPU
 * comment in RF_solvent_interaction_innerloop verbatim), and no self-term
 * (its distance-independent part is left out deliberately: it sums to
 * zero and is dropped, not computed-then-cancelled).
 */
template <math::boundary_enum BOUNDARY>
__global__ void rf_excluded_solvent_kernel(
    const int* __restrict__ chargegroup,
    unsigned num_solute_chargegroups,
    unsigned num_chargegroups,
    math::CuVArray::View pos,
    const FPL_TYPE* __restrict__ charge,
    const unsigned* __restrict__ atom_energy_group,
    gpu::NbSimParams nb,
    gpu::Periodicity<BOUNDARY> periodicity,
    double* e_crf_total) {

    const unsigned num_groups = nb.num_energy_groups;

    for (unsigned cg = num_solute_chargegroups + blockIdx.x * blockDim.x + threadIdx.x;
         cg < num_chargegroups; cg += blockDim.x * gridDim.x) {

        const unsigned a_from = static_cast<unsigned>(chargegroup[cg]);
        const unsigned a_to   = static_cast<unsigned>(chargegroup[cg + 1]);

        for (unsigned a1 = a_from; a1 < a_to; ++a1) {
            for (unsigned a2 = a1 + 1; a2 < a_to; ++a2) {
                const FPL3_TYPE rvec = periodicity.nearest_image(pos(a1), pos(a2));
                const FPL_TYPE dist2 = abs2(rvec);
                const FPL_TYPE q = charge[a1] * charge[a2];

                // No crf_cut term here -- see this file's/header's doc
                // comment, matches RF_solvent_interaction_innerloop
                // exactly (nonbonded_innerloop.cc).
                const FPL_TYPE e_crf = -q * nb.four_pi_eps_i * nb.crf_2cut3i * dist2;

                const unsigned eg1 = atom_energy_group[a1];
                const unsigned eg2 = atom_energy_group[a2];
                atomicAdd(&e_crf_total[eg1 * num_groups + eg2], static_cast<double>(e_crf));
            }
        }
    }
}

} // namespace gpu

void gpu::launch_rf_excluded(
    const int* rf_excl_ptr,
    const int* rf_excl_list,
    unsigned num_solute_atoms,
    const int* chargegroup,
    unsigned num_solute_chargegroups,
    unsigned num_chargegroups,
    math::CuVArray::View pos,
    const FPL_TYPE* charge,
    const unsigned* atom_energy_group,
    gpu::NbSimParams nb,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* force,
    double* e_crf_total,
    double* virial_total,
    cudaStream_t stream) {

    constexpr unsigned kThreads = 128;

    const unsigned num_solute_blocks =
        (num_solute_atoms + kThreads - 1) / kThreads;

    const unsigned num_solvent_chargegroups = num_chargegroups - num_solute_chargegroups;
    const unsigned num_solvent_blocks =
        (num_solvent_chargegroups + kThreads - 1) / kThreads;

    switch (boundary) {
        case math::vacuum: {
            const gpu::Periodicity<math::vacuum> periodicity(box);
            if (num_solute_atoms > 0) {
                gpu::rf_excluded_solute_kernel<math::vacuum><<<num_solute_blocks, kThreads, 0, stream>>>(
                    rf_excl_ptr, rf_excl_list, num_solute_atoms,
                    pos, charge, atom_energy_group, nb, periodicity,
                    force, e_crf_total, virial_total);
            }
            if (num_solvent_chargegroups > 0) {
                gpu::rf_excluded_solvent_kernel<math::vacuum><<<num_solvent_blocks, kThreads, 0, stream>>>(
                    chargegroup, num_solute_chargegroups, num_chargegroups,
                    pos, charge, atom_energy_group, nb, periodicity, e_crf_total);
            }
            break;
        }
        case math::rectangular: {
            const gpu::Periodicity<math::rectangular> periodicity(box);
            if (num_solute_atoms > 0) {
                gpu::rf_excluded_solute_kernel<math::rectangular><<<num_solute_blocks, kThreads, 0, stream>>>(
                    rf_excl_ptr, rf_excl_list, num_solute_atoms,
                    pos, charge, atom_energy_group, nb, periodicity,
                    force, e_crf_total, virial_total);
            }
            if (num_solvent_chargegroups > 0) {
                gpu::rf_excluded_solvent_kernel<math::rectangular><<<num_solvent_blocks, kThreads, 0, stream>>>(
                    chargegroup, num_solute_chargegroups, num_chargegroups,
                    pos, charge, atom_energy_group, nb, periodicity, e_crf_total);
            }
            break;
        }
        case math::triclinic: {
            const gpu::Periodicity<math::triclinic> periodicity(box);
            if (num_solute_atoms > 0) {
                gpu::rf_excluded_solute_kernel<math::triclinic><<<num_solute_blocks, kThreads, 0, stream>>>(
                    rf_excl_ptr, rf_excl_list, num_solute_atoms,
                    pos, charge, atom_energy_group, nb, periodicity,
                    force, e_crf_total, virial_total);
            }
            if (num_solvent_chargegroups > 0) {
                gpu::rf_excluded_solvent_kernel<math::triclinic><<<num_solvent_blocks, kThreads, 0, stream>>>(
                    chargegroup, num_solute_chargegroups, num_chargegroups,
                    pos, charge, atom_energy_group, nb, periodicity, e_crf_total);
            }
            break;
        }
        default:
            // Unsupported boundary -- matches CUDA_Pairlist_Algorithm::init()'s
            // vacuum/rectangular-only v1 scope; no kernel launch.
            break;
    }
    CUDA_CHECK_ERROR("rf_excluded_kernel");
}
