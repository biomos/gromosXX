#pragma once

#include <cuda_runtime.h>
#include <vector>
#include <stdexcept>
#include <string>
#include <unordered_map>


namespace gpu {
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
         *
         * Ensures that all created streams are properly destroyed.
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
         * @throws std::invalid_argument if the device ID is invalid.
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
         * @throws std::runtime_error if stream creation fails.
         */
        cudaStream_t create_stream();

        /**
         * @brief Destroy a CUDA stream.
         * @param stream The CUDA stream to destroy.
         * @throws std::runtime_error if stream destruction fails.
         */
        void destroy_stream(cudaStream_t stream);

        /**
         * @brief Get a human-readable description of a CUDA device.
         * @param device_id The ID of the device.
         * @return A string describing the device.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        std::string get_device_description(int device_id) const;

        /**
         * @brief Automatically select the best CUDA device based on properties.
         * @return The ID of the selected device.
         * @throws std::runtime_error if no suitable device is found.
         */
        int select_best_device() const;

        /**
         * @brief Reset the active device to its default state.
         * @param device_id The ID of the device to reset.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        void reset_device(int device_id);

    private:
        /**
         * @brief Validate a device ID.
         * @param device_id The ID of the device to validate.
         * @throws std::invalid_argument if the device ID is invalid.
         */
        void validate_device_id(int device_id) const;

        int device_count_; ///< The number of available CUDA devices.
        int active_device_id_; ///< The ID of the currently active CUDA device.
        std::unordered_map<cudaStream_t, int> stream_map_; ///< Tracks streams and their associated devices.
    };
}