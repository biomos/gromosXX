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
 * @file constraints.h
 * constraints algorithms
 */

#ifndef INCLUDED_CUKERNEL_CONSTRAINTS_H
#define INCLUDED_CUKERNEL_CONSTRAINTS_H
#include "parameter.h"

namespace cuda {
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
        cukernel::simulation_parameter * dev_params,
        int *shake_fail_mol
);

}

#endif
