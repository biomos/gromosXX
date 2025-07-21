

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include "gpu/cuda/cuheader.h"
#include "utils.h"

void gpu::check_cuda_call_error(cudaError_t result, const char* file, int line, const std::string& call) {
    if (result != cudaSuccess) {
        std::ostringstream oss;
        oss << "CUDA Error: " << cudaGetErrorString(result) << "\n"
            << "  at " << file << ":" << line << "\n"
            << "  during: " << call;
        io::messages.add(oss.str(), "gpu", io::message::error);
    }
}

void gpu::check_cuda_last_error(const char* err_msg, const char* file, int line) {
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        std::ostringstream oss;
        oss << "CUDA Error detected at " << file << ":" << line << " - " << err_msg
                << ": " << cudaGetErrorString(error) << std::endl;
        io::messages.add(oss.str(), "gpu", io::message::error);
    }
}

void gpu::print_device_info() {
    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));

    if (device_count == 0) {
        std::cout << "No CUDA devices available.\n";
        return;
    }

    std::cout << "CUDA Devices:\n";
    for (int i = 0; i < device_count; ++i) {
        cudaDeviceProp properties;
        CUDA_CHECK(cudaGetDeviceProperties(&properties, i));

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

std::string gpu::get_active_device_name() {
    int device_id = 0;
    CUDA_CHECK(cudaGetDevice(&device_id));

    cudaDeviceProp properties;
    CUDA_CHECK(cudaGetDeviceProperties(&properties, device_id));

    return std::string(properties.name);
}