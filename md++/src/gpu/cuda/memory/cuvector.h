/**
 * @file cuvector.h
 * @author poliak
 * CUDA-capable light-weight versions of std::vector
 */

#pragma once

#include <cuda_runtime.h>
#include <stdexcept>
#include <utility>
#include <new>
#include <iostream>

#include "cumallocator.h"
#include "cuhallocator.h"
#include "math/gmath.h"

namespace gpu
{
    /**
     * @brief Generic GPU view of vectors and arrays
     * 
     */
    template <typename ContainerT>
    class DeviceView {
        public:
            using value_type = typename ContainerT::value_type;

            __host__ DeviceView(value_type* ptr, size_t sz) 
                : m_data(ptr), m_size(sz) {}

            HOSTDEVICE value_type& operator[](size_t i) const {
                assert(i < m_size);
                return m_data[i];
            }
            HOSTDEVICE const value_type& operator[](size_t i) const {
                assert(i < m_size);
                return data[i];
            }

            class iterator {
                public:
                    HOSTDEVICE iterator(value_type* p = nullptr) : m_ptr(p) {}
                    HOSTDEVICE value_type& operator*() const { return *m_ptr; }
                    HOSTDEVICE iterator& operator++() { ++m_ptr; return *this; }
                    HOSTDEVICE bool operator==(const iterator& o) const { return m_ptr == o.m_ptr; }
                    HOSTDEVICE bool operator!=(const iterator& o) const { return m_ptr != o.m_ptr; }
                private:
                    value_type* m_ptr = nullptr;
            };

            class const_iterator {
                public:
                    HOSTDEVICE const_iterator(const value_type* p = nullptr) : m_ptr(p) {}
                    HOSTDEVICE value_type& operator*() const { return *m_ptr; }
                    HOSTDEVICE iterator& operator++() { ++m_ptr; return *this; }
                    HOSTDEVICE bool operator==(const iterator& o) const { return m_ptr == o.m_ptr; }
                    HOSTDEVICE bool operator!=(const iterator& o) const { return m_ptr != o.m_ptr; }
                private:
                    value_type* m_ptr = nullptr;
            };

            __host__ __device__ iterator begin() const { return iterator(m_data); }
            __host__ __device__ iterator end()   const { return iterator(m_data + m_size); }
            __host__ __device__ const_iterator begin() const { return const_iterator(m_data); }
            __host__ __device__ const_iterator end()   const { return const_iterator(m_data + m_size); }
            __host__ __device__ bool     empty() const { return m_size == 0; }
        private:
            
            value_type* m_data = nullptr;
            size_t m_size = 0;
    };
    template <template <typename> class ContainerT, typename T>
    struct UnifiedMemoryArray : public ContainerT<T, CuMAllocator<T>> {
        using Base = ContainerT<T, CuMAllocator<T>>;
        using View = DeviceView<Base>;

        __host__ View view() const { return View{this->data(), this->size()}; }
    };

    /**
     * @brief CUDA-accessible variant of std::vector
     * 
     * @tparam T 
     */
    template <typename T>
    using cuvector = UnifiedMemoryArray<std::vector, T>;

    /**
     * @brief Pinned-memory variant of std::vector
     * 
     * @tparam T 
     */
    template <typename T>
    using cuhvector = std::vector< T, CuHAllocator<T> >;
}

namespace math
{
    /**
     * @brief CUDA-accessible variant of math::VArrayT
     * 
     * @tparam T 
     */
    template <typename T>
    using CuVArrayT = UnifiedMemoryArray<math::VArrayT, T>;

    /**
     * @brief CUDA-accessible variant of math::VArray
     * 
     * @tparam T 
     */
    using CuVArray = CuVArrayT<FPL3_TYPE>;
}
