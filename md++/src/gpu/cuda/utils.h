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
 * @file utils.h
 * @brief A collection of utility functions for GPU-related operations.
 */

#pragma once

#include <cuda_runtime.h>
#include <string>

#define CHECK(call) gpu::check_cuda_error((call), __FILE__, __LINE__, #call)

namespace gpu {
    /**
     * @brief Check for CUDA errors and throw an exception if an error is detected.
     * @param err The CUDA error code to check.
     * @param context A string describing the context of the error (e.g., function name).
     * @throws std::runtime_error if a CUDA error is detected.
     */
    void check_cuda_error(cudaError_t err, const char* file, int line, const std::string& call);

    /**
     * @brief Query and print information about all available CUDA devices.
     */
    void print_device_info();

    /**
     * @brief Synchronize the device and measure the elapsed time of a GPU operation.
     * @param start The CUDA event marking the start of the operation.
     * @param stop The CUDA event marking the end of the operation.
     * @return The elapsed time in milliseconds.
     */
    inline float measure_elapsed_time(cudaEvent_t start, cudaEvent_t stop) {
        CHECK(cudaEventSynchronize(stop));
        float milliseconds = 0.0f;
        CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
        return milliseconds;
    }

    /**
     * @brief Create and return a CUDA event.
     * @return A CUDA event.
     */
    inline cudaEvent_t create_event() {
        cudaEvent_t event;
        CHECK(cudaEventCreate(&event));
        return event;
    };

    /**
     * @brief Destroy a CUDA event.
     * @param event The CUDA event to destroy.
     */
    inline void destroy_event(cudaEvent_t event) {
        CHECK(cudaEventDestroy(event));
    };

    /**
     * @brief Get the name of the currently active CUDA device.
     * @return The name of the active CUDA device.
     */
    std::string get_active_device_name();

} // namespace gpu
