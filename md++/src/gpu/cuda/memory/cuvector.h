/**
 * @file cuvector.h
 * @author poliak
 * CUDA-capable light-weight versions of std::vector
 */

#pragma once

#include "cumallocator.h"
#include "cuhallocator.h"
#include "math/gmath.h"

namespace gpu
{
    // template <typename T>
    // using cuvector = std::vector< T, CuMAllocator<T> >;
    template <typename T>
    using cuhvector = std::vector< T, CuHAllocator<T> >;

    template <typename T>
    class cuvector {
    public:
        using value_type = T;
        using allocator_type = CuMAllocator<T>;
        using size_type = std::size_t;

        // Device-accessible iterator
        struct iterator {
            T* ptr;

            __host__ __device__
            iterator(T* p = nullptr) : ptr(p) {}

            __host__ __device__
            T& operator*() const { return *ptr; }

            __host__ __device__
            iterator& operator++() { ++ptr; return *this; }

            __host__ __device__
            iterator operator+(size_type n) const { return iterator(ptr + n); }

            __host__ __device__
            bool operator!=(const iterator& other) const { return ptr != other.ptr; }

            __host__ __device__
            T* operator->() const { return ptr; }
        };

        // Device-accessible const iterator
        struct const_iterator {
            const T* ptr;

            __host__ __device__
            const_iterator(const T* p = nullptr) : ptr(p) {}

            __host__ __device__
            const T& operator*() const { return *ptr; }

            __host__ __device__
            const_iterator& operator++() { ++ptr; return *this; }

            __host__ __device__
            const_iterator operator+(size_type n) const { return const_iterator(ptr + n); }

            __host__ __device__
            bool operator!=(const const_iterator& other) const { return ptr != other.ptr; }

            __host__ __device__
            const T* operator->() const { return ptr; }
        };

        // Constructors / destructors
        __host__ cuvector() : m_data(nullptr), m_size(0), m_capacity(0) {}
        
        __host__ explicit cuvector(size_type n)
            : m_size(n), m_capacity(n) 
        {
            m_data = m_allocator.allocate(n);
            for (size_type i = 0; i < n; ++i) {
                m_allocator.construct(&m_data[i]);
            }
        }

        __host__ ~cuvector() {
            for (size_type i = 0; i < m_size; ++i)
                m_allocator.destroy(&m_data[i]);
            m_allocator.deallocate(m_data, m_capacity);
        }

        // Size
        __host__ __device__
        size_type size() const { return m_size; }

        // Access
        __host__ __device__
        T& operator[](size_type i) { return m_data[i]; }

        __host__ __device__
        const T& operator[](size_type i) const { return m_data[i]; }

        __host__ __device__
        T* data() { return m_data; }

        __host__ __device__
        const T* data() const { return m_data; }

        // Iterators
        __host__ __device__
        iterator begin() { return iterator(m_data); }

        __host__ __device__
        iterator end() { return iterator(m_data + m_size); }

        __host__ __device__
        const_iterator begin() const { return const_iterator(m_data); }

        __host__ __device__
        const_iterator end() const { return const_iterator(m_data + m_size); }

        __host__ __device__
        const_iterator cbegin() const { return const_iterator(m_data); }

        __host__ __device__
        const_iterator cend() const { return const_iterator(m_data + m_size); }

        // Capacity
        __host__
        void reserve(size_type new_cap) {
            if (new_cap <= m_capacity) return;

            T* new_data = m_allocator.allocate(new_cap);
            for (size_type i = 0; i < m_size; ++i) {
                m_allocator.construct(&new_data[i], std::move(m_data[i]));
                m_allocator.destroy(&m_data[i]);
            }
            m_allocator.deallocate(m_data, m_capacity);
            m_data = new_data;
            m_capacity = new_cap;
        }

        __host__
        void resize(size_type new_size) {
            if (new_size > m_capacity) reserve(new_size);
            for (size_type i = m_size; i < new_size; ++i)
                m_allocator.construct(&m_data[i]);
            for (size_type i = new_size; i < m_size; ++i)
                m_allocator.destroy(&m_data[i]);
            m_size = new_size;
        }

        __host__
        void push_back(const T& value) {
            if (m_size >= m_capacity) {
                size_type new_cap = (m_capacity == 0) ? 1 : m_capacity * 2;
                reserve(new_cap);
            }
            m_allocator.construct(&m_data[m_size], value);
            ++m_size;
        }

        __host__
        void push_back(T&& value) {
            if (m_size >= m_capacity) {
                size_type new_cap = (m_capacity == 0) ? 1 : m_capacity * 2;
                reserve(new_cap);
            }
            m_allocator.construct(&m_data[m_size], std::move(value));
            ++m_size;
        }

        __host__
        void clear() noexcept {
            for (size_type i = 0; i < m_size; ++i) {
                std::destroy_at(m_data + i);
            }
            m_size = 0;
        }

    private:
        allocator_type m_allocator;
        T* m_data;
        size_type m_size;
        size_type m_capacity;
    };
}

namespace math
{
    /**
     * CUDA allocated managed VArray (double precision)
     */
    // typedef VArrayT<Vec, gpu::CuMAllocator > CuDVArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    // typedef VArrayT<Vecf, gpu::CuMAllocator > CuVArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    // typedef VArrayT<float3, gpu::CuMAllocator > CuVArray;
    typedef gpu::cuvector<float3> CuVArray;

    /**
     * CUDA allocated managed SArray (single precision)
     */
    typedef gpu::cuvector<float> CuSArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    // typedef VArrayT<double3, gpu::CuMAllocator > CuD3VArray;
}
