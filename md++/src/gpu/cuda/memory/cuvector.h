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

            HOSTDEVICE DeviceView(value_type* ptr, size_t sz) 
                : m_data(ptr), m_size(sz) {}

            HOSTDEVICE value_type& operator[](size_t i) {
                return m_data[i];
            }
            HOSTDEVICE const value_type& operator[](size_t i) const {
                return m_data[i];
            }

            HOSTDEVICE value_type& operator()(size_t i) {
                assert(i < m_size);
                return this->operator[](i);
            }
            HOSTDEVICE const value_type& operator()(size_t i) const {
                assert(i < m_size);
                return this->operator[](i);
            }

            HOSTDEVICE value_type& at(size_t i) {
                if (i >= m_size) {
                    #ifndef __CUDA_ARCH__   // host side
                        throw std::out_of_range("DeviceView::at() index out of range");
                    #else                   // device side
                        printf("DeviceView::at() index out of range: %zu (size=%zu)\n", i, m_size);
                        asm("trap;"); // abort thread
                    #endif
                }
                return m_data[i];
            }

            HOSTDEVICE bool empty() const { return m_size == 0; }
            HOSTDEVICE size_t size() const { return m_size; }

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
                    HOSTDEVICE const value_type& operator*() const { return *m_ptr; }
                    HOSTDEVICE const_iterator& operator++() { ++m_ptr; return *this; }
                    HOSTDEVICE bool operator==(const const_iterator& o) const { return m_ptr == o.m_ptr; }
                    HOSTDEVICE bool operator!=(const const_iterator& o) const { return m_ptr != o.m_ptr; }
                private:
                    const value_type* m_ptr = nullptr;
            };

            // Non-const begin/end
            HOSTDEVICE iterator begin() { return iterator(m_data); }
            HOSTDEVICE iterator end()   { return iterator(m_data + m_size); }

            // Const begin/end
            HOSTDEVICE const_iterator begin() const { return const_iterator(m_data); }
            HOSTDEVICE const_iterator end()   const { return const_iterator(m_data + m_size); }

        private:
            value_type* const m_data;
            const size_t m_size;
    };

    /**
     * @brief CUDA-accessible variant of std::vector and math::VArray
     * 
     * @tparam ContainerT 
     * @tparam T 
     */
    template <template <typename, typename> class ContainerT, typename T>
    struct UnifiedMemoryArray : public ContainerT<T, CuMAllocator<T>> {
        using Base = ContainerT<T, CuMAllocator<T>>;
        using View = DeviceView<Base>;

        __host__ View view() { return View{this->data(), this->size()}; }
        __host__ const View const_view() const { return View{this->data(), this->size()}; }
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
    using CuVArrayT = gpu::UnifiedMemoryArray<math::VArrayT, T>;

    /**
     * @brief CUDA-accessible variant of math::VArray
     * 
     * @tparam T 
     */
    using CuVArray = CuVArrayT<FPL3_TYPE>;
}
