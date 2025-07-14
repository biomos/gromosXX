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
 * @file constraints.cu
 * contains the GPU methods of constraints.
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "gpu.h"
#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/manager/cuda_manager.tcc"

#include "../../../io/print_block.h"

#include "constraints.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE algorithm

double gpu::cuda::algorithm::remove_com_translation(CudaManager& cuda_manager,
                                const topology::Topology& topo,
                                const configuration::Configuration& conf,
                                const simulation::Simulation& sim,
                                bool remove_trans) {
    // CUDA-specific implementation
    cuda_manager.synchronize_all();
    auto pos = cuda_manager.create_cuvector<float3>(0, topo.num_atoms());
    return 0;
}

double gpu::cuda::algorithm::remove_com_rotation(CudaManager& cuda_manager,
                                const topology::Topology& topo,
                                const configuration::Configuration& conf,
                                const simulation::Simulation& sim,
                                bool remove_trans) {
    // CUDA-specific implementation
    cuda_manager.synchronize_all();
    auto pos = cuda_manager.create_cuvector<float3>(0, topo.num_atoms());
    return 0;
}