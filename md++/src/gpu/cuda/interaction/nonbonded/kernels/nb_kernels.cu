/**
 * @file nb_kernels.cu
 * GPU kernels for vanilla LJ + reaction-field (CRF) nonbonded forces.
 */

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/float3.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/utils.h"
#include "math/boundary.h"

#include "nb_kernels.h"

// ─────────────────────────────────────────────────────────────────
// Kernel: zero force array
// ─────────────────────────────────────────────────────────────────

__global__ void gpu::zero_forces_kernel(FPL3_TYPE* force, unsigned num_atoms) {
    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (unsigned i = idx; i < num_atoms; i += blockDim.x * gridDim.x)
        force[i] = make_FPL3(FPL_TYPE(0));
}

// ─────────────────────────────────────────────────────────────────
// Kernel: LJ-CRF pairwise forces
// ─────────────────────────────────────────────────────────────────

template <math::boundary_enum BOUNDARY>
__global__ void gpu::lj_crf_forces_kernel(
    const FPL3_TYPE* __restrict__ pos,
    FPL3_TYPE*                    force,
    const int*       __restrict__ iac,
    const FPL_TYPE*  __restrict__ charge,
    const uint2*     __restrict__ pairs,
    unsigned                      num_pairs,
    gpu::LJParamView              lj,
    gpu::NbSimParams              nb,
    gpu::Periodicity<BOUNDARY>    periodicity)
{
    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    for (unsigned k = idx; k < num_pairs; k += blockDim.x * gridDim.x) {
        const unsigned i = pairs[k].x;
        const unsigned j = pairs[k].y;

        // r = nearest_image(pos[i], pos[j])  → vector j→i (GROMOS convention)
        FPL3_TYPE r = periodicity.nearest_image(pos[i], pos[j]);
        const FPL_TYPE dist2 = r.x*r.x + r.y*r.y + r.z*r.z;

        // ── Lennard-Jones ────────────────────────────────────────
        const FPL_TYPE dist6i = FPL_TYPE(1) / (dist2 * dist2 * dist2);
        const FPL2_TYPE lj_p  = lj.get(iac[i], iac[j]);
        const FPL_TYPE c6     = lj_p.x;
        const FPL_TYPE c12    = lj_p.y;

        // f_lj = (12*c12/r^14 - 6*c6/r^8)  – positive = repulsive in j→i dir
        FPL_TYPE f = (FPL_TYPE(12) * c12 * dist6i - FPL_TYPE(6) * c6) * dist6i / dist2;

        // ── Reaction-field Coulomb ────────────────────────────────
        const FPL_TYPE ri  = rsqrtf(dist2);          // 1/r  (float rsqrt is faster)
        const FPL_TYPE ri3 = ri / dist2;             // 1/r^3
        const FPL_TYPE q   = charge[i] * charge[j]; // charge product

        f += q * nb.four_pi_eps_i * (ri3 + FPL_TYPE(2) * nb.crf_2cut3i);

        // ── Force accumulation (Newton 3rd law) ──────────────────
        // F_i += f * r,  F_j -= f * r
        atomicAdd(&force[i].x,  f * r.x);
        atomicAdd(&force[i].y,  f * r.y);
        atomicAdd(&force[i].z,  f * r.z);
        atomicAdd(&force[j].x, -f * r.x);
        atomicAdd(&force[j].y, -f * r.y);
        atomicAdd(&force[j].z, -f * r.z);
    }
}

// Explicit template instantiations (required for separate compilation)
template __global__ void gpu::lj_crf_forces_kernel<math::vacuum>(
    const FPL3_TYPE*, FPL3_TYPE*, const int*, const FPL_TYPE*,
    const uint2*, unsigned, gpu::LJParamView, gpu::NbSimParams,
    gpu::Periodicity<math::vacuum>);

template __global__ void gpu::lj_crf_forces_kernel<math::rectangular>(
    const FPL3_TYPE*, FPL3_TYPE*, const int*, const FPL_TYPE*,
    const uint2*, unsigned, gpu::LJParamView, gpu::NbSimParams,
    gpu::Periodicity<math::rectangular>);

template __global__ void gpu::lj_crf_forces_kernel<math::triclinic>(
    const FPL3_TYPE*, FPL3_TYPE*, const int*, const FPL_TYPE*,
    const uint2*, unsigned, gpu::LJParamView, gpu::NbSimParams,
    gpu::Periodicity<math::triclinic>);

// ─────────────────────────────────────────────────────────────────
// Host-side launch wrappers
// ─────────────────────────────────────────────────────────────────

static constexpr unsigned BLOCK = 256;

void gpu::launch_zero_forces(FPL3_TYPE* d_force, unsigned num_atoms,
                              cudaStream_t stream)
{
    if (num_atoms == 0) return;
    const unsigned grid = (num_atoms + BLOCK - 1) / BLOCK;
    zero_forces_kernel<<<grid, BLOCK, 0, stream>>>(d_force, num_atoms);
    CUDA_CHECK_ERROR("zero_forces_kernel");
}

void gpu::launch_lj_crf_forces(
    const FPL3_TYPE* pos,
    FPL3_TYPE*       force,
    const int*       iac,
    const FPL_TYPE*  charge,
    const uint2*     pairs,
    unsigned         num_pairs,
    gpu::LJParamView lj,
    gpu::NbSimParams nb,
    math::boundary_enum boundary,
    math::Box        box,
    cudaStream_t     stream)
{
    if (num_pairs == 0) return;
    const unsigned grid = (num_pairs + BLOCK - 1) / BLOCK;

    switch (boundary) {
        case math::vacuum:
            lj_crf_forces_kernel<math::vacuum>
                <<<grid, BLOCK, 0, stream>>>(
                    pos, force, iac, charge, pairs, num_pairs, lj, nb,
                    gpu::Periodicity<math::vacuum>(box));
            break;
        case math::rectangular:
            lj_crf_forces_kernel<math::rectangular>
                <<<grid, BLOCK, 0, stream>>>(
                    pos, force, iac, charge, pairs, num_pairs, lj, nb,
                    gpu::Periodicity<math::rectangular>(box));
            break;
        case math::triclinic:
            lj_crf_forces_kernel<math::triclinic>
                <<<grid, BLOCK, 0, stream>>>(
                    pos, force, iac, charge, pairs, num_pairs, lj, nb,
                    gpu::Periodicity<math::triclinic>(box));
            break;
        default:
            // unsupported boundary – fall through without kernel launch
            break;
    }
    CUDA_CHECK_ERROR("lj_crf_forces_kernel");
}
