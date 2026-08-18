/**
 * @file cuda_lj_params.cu
 * Builds the flat GPU LJ parameter matrix from the CPU nonbonded parameter table.
 */

#include "gpu/cuda/cuheader.h"

#include "stdheader.h"
#include "interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "cuda_lj_params.h"

void gpu::LJParams::init(interaction::Nonbonded_Parameter& params) {
    // Nonbonded_Parameter::lj_parameter() (no-arg, full-matrix accessor)
    // isn't const, hence the non-const reference here.
    const auto& matrix = params.lj_parameter();
    num_types = static_cast<unsigned>(matrix.size());

    const unsigned total = num_types * num_types;
    c6.resize(total);
    c12.resize(total);
    cs6.resize(total);
    cs12.resize(total);

    for (unsigned i = 0; i < num_types; ++i) {
        for (unsigned j = 0; j < num_types; ++j) {
            const unsigned idx = i * num_types + j;
            c6  [idx] = static_cast<FPL_TYPE>(matrix[i][j].c6);
            c12 [idx] = static_cast<FPL_TYPE>(matrix[i][j].c12);
            cs6 [idx] = static_cast<FPL_TYPE>(matrix[i][j].cs6);
            cs12[idx] = static_cast<FPL_TYPE>(matrix[i][j].cs12);
        }
    }
}
