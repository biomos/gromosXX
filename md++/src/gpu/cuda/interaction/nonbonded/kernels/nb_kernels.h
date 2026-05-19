/**
 * @file nb_kernels.h
 * GPU kernel declarations for vanilla LJ-CRF nonbonded interactions.
 *
 * Two kernels are provided:
 *   - zero_forces_kernel  : zero the force array before accumulation
 *   - lj_crf_forces_kernel: one thread per pair, LJ + reaction-field forces
 *
 * Host-side launch wrappers handle boundary dispatch.
 */

#pragma once

#include "gpu/cuda/cuhostdevice.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "math/boundary.h"

namespace gpu {

/** Zero the GPU force array (one thread per atom). */
__global__ void zero_forces_kernel(FPL3_TYPE* force, unsigned num_atoms);

/**
 * @brief Compute LJ + reaction-field forces for a flat atom-pair list.
 *
 * One CUDA thread processes one pair (i, j).
 * Forces are accumulated with atomicAdd (no reduction needed).
 *
 * Convention (matching GROMOS CPU):
 *   r   = nearest_image(pos[i], pos[j])  (vector j→i)
 *   F_i += f_scalar * r
 *   F_j -= f_scalar * r
 *
 * @tparam BOUNDARY  Periodic boundary type (vacuum / rectangular / triclinic)
 */
template <math::boundary_enum BOUNDARY>
__global__ void lj_crf_forces_kernel(
    const FPL3_TYPE* __restrict__ pos,
    FPL3_TYPE*                    force,
    const int*       __restrict__ iac,
    const FPL_TYPE*  __restrict__ charge,
    const uint2*     __restrict__ pairs,
    unsigned                      num_pairs,
    gpu::LJParamView              lj,
    gpu::NbSimParams              nb,
    gpu::Periodicity<BOUNDARY>    periodicity);

// -----------------------------------------------------------------
// Host-side launch wrappers (dispatch on boundary type at runtime)
// -----------------------------------------------------------------

/** Zero the force array on device. */
void launch_zero_forces(FPL3_TYPE* d_force, unsigned num_atoms,
                        cudaStream_t stream = nullptr);

/**
 * @brief Launch lj_crf_forces_kernel for the given boundary type.
 *
 * @param pos          Device pointer to positions
 * @param force        Device pointer to forces (accumulated in-place)
 * @param iac          Device pointer to integer atom codes
 * @param charge       Device pointer to charges
 * @param pairs        Device pointer to flat pair list (uint2 = {i, j})
 * @param num_pairs    Number of pairs
 * @param lj           LJ parameter view
 * @param nb           Simulation nonbonded parameters
 * @param boundary     Boundary type at runtime
 * @param box          Host-side box (used to build Periodicity object)
 * @param stream       CUDA stream (default = default stream)
 */
void launch_lj_crf_forces(
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
    cudaStream_t     stream = nullptr);

} // namespace gpu
