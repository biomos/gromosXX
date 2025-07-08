#pragma once

#include <cuda_runtime.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

/**
 * @namespace gpu
 * @brief A collection of utility functions for GPU-related operations.
 */

#define CHECK(call) gpu::check_cuda_error((call), __FILE__, __LINE__, #call)

namespace gpu {
    /**
     * @brief Check for CUDA errors and throw an exception if an error is detected.
     * @param err The CUDA error code to check.
     * @param context A string describing the context of the error (e.g., function name).
     * @throws std::runtime_error if a CUDA error is detected.
     */
    inline void check_cuda_error(cudaError_t err, const char* file, int line, const std::string& call) {
        if (err != cudaSuccess) {
            std::ostringstream oss;
            oss << "CUDA Error: " << cudaGetErrorString(err) << "\n"
                << "  at " << file << ":" << line << "\n"
                << "  during: " << call;
            throw std::runtime_error(oss.str());
        }
    }

    /**
     * @brief Query and print information about all available CUDA devices.
     */
    inline void print_device_info() {
        int device_count = 0;
        CHECK(cudaGetDeviceCount(&device_count));

        if (device_count == 0) {
            std::cout << "No CUDA devices available.\n";
            return;
        }

        std::cout << "CUDA Devices:\n";
        for (int i = 0; i < device_count; ++i) {
            cudaDeviceProp properties;
            CHECK(cudaGetDeviceProperties(&properties, i));

            std::cout << "Device " << i << ": " << properties.name << "\n"
                    << "  Compute Capability: " << properties.major << "." << properties.minor << "\n"
                    << "  Total Global Memory: " << (properties.totalGlobalMem / (1024 * 1024)) << " MB\n"
                    << "  Multiprocessors: " << properties.multiProcessorCount << "\n"
                    << "  Max Threads per Block: " << properties.maxThreadsPerBlock << "\n"
                    << "  Max Threads per SM: " << properties.maxThreadsPerMultiProcessor << "\n"
                    << "  Warp Size: " << properties.warpSize << "\n"
                    << "  Max Grid Size: (" << properties.maxGridSize[0] << ", "
                    << properties.maxGridSize[1] << ", " << properties.maxGridSize[2] << ")\n"
                    << "  Max Threads Dim: (" << properties.maxThreadsDim[0] << ", "
                    << properties.maxThreadsDim[1] << ", " << properties.maxThreadsDim[2] << ")\n"
                    << std::endl;
        }
    }

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
    }

    /**
     * @brief Destroy a CUDA event.
     * @param event The CUDA event to destroy.
     */
    inline void destroy_event(cudaEvent_t event) {
        CHECK(cudaEventDestroy(event));
    }

    /**
     * @brief Get the name of the currently active CUDA device.
     * @return The name of the active CUDA device.
     */
    inline std::string get_active_device_name() {
        int device_id = 0;
        CHECK(cudaGetDevice(&device_id));

        cudaDeviceProp properties;
        CHECK(cudaGetDeviceProperties(&properties, device_id));

        return std::string(properties.name);
    }
} // namespace gpu