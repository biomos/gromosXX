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
 * @file configuration_struct.cu
 * Implementation of the light-weight configuration struct for GPU
 */

#include "stdheader.h"

#include "configuration/configuration_global.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "configuration/configuration.h"
#include "configuration/mesh.h"
#include "configuration/influence_function.h"
#include "simulation/simulation.h"
#include "simulation/multibath.h"
#include "simulation/parameter.h"

#include "math/periodicity.h"
#include "math/boundary_checks.h"
#include "util/template_split.h"

#include "configuration_struct.h"

#include "gpu/cuda/utils.h"

void gpu::Configuration::copy_to_device(configuration::Configuration& conf) {
    CUDA_CHECK_ERROR("At gpu::Configuration::copy_to_device");
    const size_t num_atoms = conf.current().pos.size();
    using Vec = typename decltype(conf.current().pos)::value_type;

    static_assert(std::is_convertible<Vec, FPL3_TYPE>::value,
                  "Vec must be convertible to FPL3_TYPE");

    auto convert_and_copy = [num_atoms](const auto& src, auto& dst) {
        CUDA_CHECK_ERROR("Before resize");
        dst.resize(num_atoms);  // Resize already allocates GPU memory
        CUDA_CHECK_ERROR("After resize");
        for (size_t i = 0; i < num_atoms; ++i)
            dst[i] = static_cast<FPL3_TYPE>(src[i]);  // Copy/convert element-wise
    };

    // Current state
    convert_and_copy(conf.current().pos, current.pos);
    convert_and_copy(conf.current().vel, current.vel);
    convert_and_copy(conf.current().force, current.force);
    convert_and_copy(conf.current().constraint_force, current.constraint_force);

    // Old state
    convert_and_copy(conf.old().pos, old.pos);
    convert_and_copy(conf.old().vel, old.vel);
    convert_and_copy(conf.old().force, old.force);
    convert_and_copy(conf.old().constraint_force, old.constraint_force);

    
    CUDA_CHECK_ERROR("Before tensors");
    // Copy raw Box and tensor pointers
    CUDA_CHECK(cudaMemcpy(current.box, &conf.current().box, sizeof(Box), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(current.virial_tensor, &conf.current().virial_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(current.kinetic_energy_tensor, &conf.current().kinetic_energy_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(current.pressure_tensor, &conf.current().pressure_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(old.box, &conf.old().box, sizeof(Box), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.virial_tensor, &conf.old().virial_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.kinetic_energy_tensor, &conf.old().kinetic_energy_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.pressure_tensor, &conf.old().pressure_tensor, sizeof(FPL9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR("After tensors");
}