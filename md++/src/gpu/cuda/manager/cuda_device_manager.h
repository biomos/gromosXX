#pragma once

#include <cuda_runtime.h>
#include <vector>
#include <stdexcept>
#include <string>

/**
 * @class CudaDeviceManager
 * @brief Manages CUDA devices and streams.
 *
 * The CudaDeviceManager class provides functionality to query available CUDA devices,
 * set the active device, and manage CUDA streams. It abstracts away low-level CUDA
 * device management details and provides a clean interface for interacting with GPUs.
 */
class CudaDeviceManager {
public:
    /**
     * @brief Constructor for CudaDeviceManager.
     *
     * Queries the available CUDA devices and initializes the manager.
     * Throws an exception if no CUDA devices are available.
     */
    CudaDeviceManager();

    /**
     * @brief Destructor for CudaDeviceManager.
     */
    ~CudaDeviceManager();

    /**
     * @brief Get the number of available CUDA devices.
     * @return The number of available CUDA devices.
     */
    int get_device_count() const;

    /**
     * @brief Get the properties of a specific CUDA device.
     * @param device_id The ID of the device.
     * @return The properties of the specified device.
     */
    cudaDeviceProp get_device_properties(int device_id) const;

    /**
     * @brief Set the active CUDA device.
     * @param device_id The ID of the device to set as active.
     * @throws std::invalid_argument if the device ID is invalid.
     */
    void set_device(int device_id);

    /**
     * @brief Get the ID of the currently active CUDA device.
     * @return The ID of the currently active CUDA device.
     */
    int get_active_device() const;

    /**
     * @brief Create a CUDA stream on the active device.
     * @return The created CUDA stream.
     */
    cudaStream_t create_stream();

    /**
     * @brief Destroy a CUDA stream.
     * @param stream The CUDA stream to destroy.
     */
    void destroy_stream(cudaStream_t stream);

    /**
     * @brief Get a human-readable description of a CUDA device.
     * @param device_id The ID of the device.
     * @return A string describing the device.
     */
    std::string get_device_description(int device_id) const;

private:
    int device_count_; ///< The number of available CUDA devices.
    int active_device_id_; ///< The ID of the currently active CUDA device.
};
