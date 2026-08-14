/**
 * @file cuda_lj_params.h
 * GPU-resident Lennard-Jones parameter matrix.
 *
 * The matrix is stored in a flat, unified-memory array:
 *   c6 [iac_i * num_types + iac_j]
 *   c12[iac_i * num_types + iac_j]
 *
 * @note  Only the regular (non-1-4) parameters are stored here.
 */

#pragma once

#include "gpu/cuda/cuhostdevice.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"

namespace interaction {
    class Nonbonded_Parameter;
}

namespace gpu {

/**
 * @brief Lightweight device-accessible view of the LJ parameter matrix.
 * Passed by value into kernels.
 */
struct LJParamView {
    const FPL_TYPE* c6;
    const FPL_TYPE* c12;
    unsigned num_types;

    /** Get (c6, c12) for atom-type pair (ti, tj). */
    HOSTDEVICE FPL2_TYPE get(int ti, int tj) const {
        unsigned idx = (unsigned)ti * num_types + (unsigned)tj;
        return { c6[idx], c12[idx] };
    }
};

/**
 * @brief Owns unified-memory LJ parameter arrays and constructs them
 *        from the CPU-side @c interaction::Nonbonded_Parameter.
 */
struct LJParams {
    using View = LJParamView;

    gpu::cuvector<FPL_TYPE> c6;
    gpu::cuvector<FPL_TYPE> c12;
    unsigned num_types = 0;

    /**
     * Initialise / update from the CPU nonbonded parameter table.
     * Non-const: Nonbonded_Parameter::lj_parameter() (the no-arg,
     * full-matrix accessor this reads) isn't const either.
     */
    void init(interaction::Nonbonded_Parameter& params);

    View view() const {
        return View{ c6.data(), c12.data(), num_types };
    }
};

} // namespace gpu
