/**
 * @file constraints.h
 * constraints algorithms
 */

#ifndef CUKERNEL_CONSTRAINTS
#include "parameter.h"

namespace cudakernel {
/**
 * solve the constraints using the SETTLE algorithm
 * @param[inout] new_pos the new positions
 * @param[in] old_pos the old positions
 * @param[in] dev_params the simulation parameters
 * @param[out] shake_file_mol. The molecule for which the constraints algorithm failed. -1 for success.
 */
__global__ void kernel_CalcConstraints_Settle
(
        double3 * new_pos, double3 * old_pos,
        cudakernel::simulation_parameter * dev_params,
        int *shake_fail_mol
);

}

#endif
