
#pragma once

template <typename T>
class CudaArray {
private:
    T* m_device_ptr = nullptr;
    size_t m_size = 0;  // Number of elements, not bytes

public:
    CudaArray(T* device_ptr, size_t size)
        : m_device_ptr(device_ptr), m_size(size) {}

    ~CudaArray() {
        if (device_ptr) {
            cudaFree(device_ptr);
        }
    }

    T* get() const {
        return device_ptr;
    }

    size_t size() const {
        return m_size;
    }

    void copy_to_device(const T* host_data) {
        cudaMemcpy(device_ptr, host_data, m_size * sizeof(T), cudaMemcpyHostToDevice);
    }

    void copy_to_host(T* host_data) const {
        cudaMemcpy(host_data, device_ptr, m_size * sizeof(T), cudaMemcpyDeviceToHost);
    }
};

template <typename T>
class CudaManagedArray {
private:
    T* m_device_ptr = nullptr;
    size_t m_size = 0;  // Number of elements, not bytes

public:
    CudaArray(T* device_ptr, size_t size)
        : m_device_ptr(device_ptr), m_size(size) {}

    ~CudaArray() {
        if (m_device_ptr) {
            cudaFree(m_device_ptr);
        }
    }

    T* get_device_ptr() const {
        return m_device_ptr;
    }

    size_t get_size() const {
        return m_size;
    }

    void copy_to_device(const T* host_data) {
        cudaMemcpy(m_device_ptr, host_data, m_size * sizeof(T), cudaMemcpyHostToDevice);
    }

    void copy_to_host(T* host_data) const {
        cudaMemcpy(host_data, m_device_ptr, m_size * sizeof(T), cudaMemcpyDeviceToHost);
    }
};
