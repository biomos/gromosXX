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
 * @file cuheader.h
 * general header file for CUDA related code
 */

#pragma once

#include <map> // for io/message.h
#include <iostream> // for io/message.h

#include "config.h"

// #include "stdheader.h"
#include "io/message.h"

#ifdef USE_CUDA
    #include <cuda_runtime.h>
    #include "memory/cuvector.h"
    namespace gpu {
        using CUDEVPROP = cudaDeviceProp;
        using CUSTREAM = cudaStream_t;
        using CUERROR = cudaError_t;
        template <typename T>
        using CUVECTOR_T = cuvector<T>;
        template <typename T>
        using CUHVECTOR_T = cuhvector<T>;
    }
#else
    #include <vector>
    namespace gpu {
        using CUDEVPROP = void*;
        using CUSTREAM = void*;
        using CUERROR = void*;
        template <typename T>
        using CUVECTOR_T = std::vector<T>;
        template <typename T>
        using CUHVECTOR_T = std::vector<T>;
    }
#endif

namespace gpu {
    /**
     * @brief Helper function called everytime a CUDA function is called in an unsupported compilation
     * @throws io::message::critical always
     */
    template <typename T>
    auto disabled_impl(const char* file, int line, const char* func, const char* msg = nullptr)
    {
#ifdef NDEBUG
        // Release build: minimal error message
        io::messages.add(
            msg ? msg : "CUDA is disabled in this build.",
            "GPU::Cuda", io::message::critical);
#else
        // Debug build: detailed message
        io::messages.add(
            (msg ? (std::string(msg) + ". ") : "CUDA is disabled in this build. ") +
            "At " + std::string(file) + ":" + std::to_string(line) +
            " in function " + std::string(func) + "\n",
            "GPU::Cuda", io::message::critical);
#endif
        if constexpr (std::is_void_v<T>) {
            // For void, no return
            return;
        } else if constexpr (std::is_reference_v<T>) {
            // For references, return a reference to a static default value
            using RefType = std::remove_reference_t<T>;
            static RefType dummy{};
            return static_cast<T>(dummy);
        } else {
            // For normal types, return default-constructed value
            return T{};
        }
    }
}

#define DISABLED_VOID() ::gpu::disabled_impl<void>(__FILE__, __LINE__, __func__)
#define DISABLED(T) ::gpu::disabled_impl<T>(__FILE__, __LINE__, __func__)
#define DISABLED_MSG(T, msg) ::gpu::disabled_impl<T>(__FILE__, __LINE__, __func__, msg)
