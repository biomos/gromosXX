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
 * @file backend.h
 * tags for CPU and GPU variants of algorithms
 */
#pragma once

namespace gpu {
    class CudaManager;
}

namespace util {
    struct cpuBackend {};
    struct gpuBackend {};
    
    // Helper traits struct (empty by default)
    template<typename Backend>
    struct BackendData {
    // empty for all except gpu
    };

    // Specialization for gpu backend
    template<>
    struct BackendData<util::gpuBackend> {
    std::shared_ptr<::gpu::CudaManager> cuda_manager_;
    };
}