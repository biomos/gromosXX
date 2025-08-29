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
 * @file periodicity.cu
 * implementation of the periodic boundary condition functions.
 */

#include "stdheader.h"
#include "gpu/cuda/cuheader.h"

#include "math/gmath.h"
#include "math/boundary_implementation.h"

#include "algorithm/algorithm.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"

#include "topology/topology.h"

#include "configuration/configuration.h"

#include "periodicity.h"


#include "gpu/cuda/kernels/hello_world.h"
#include "gpu/cuda/kernels/periodicity.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

#define NUM_THREADS_PER_BLOCK 256
