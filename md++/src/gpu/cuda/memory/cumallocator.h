/**
 * @file Cumallocator.h
 * @author poliak
 * allocator for CUDA managed unified memory, CUDA transparent version of std::vector
 */

#pragma once

#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

#include "gpu/cuda/cuheader.h"

namespace gpu {

    template <typename T>
    class CuMAllocator {
    public:
        using value_type = T;
        using pointer = T*;
        using const_pointer = const T*;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using allocator_type = CuMAllocator<T>;

        using propagate_on_container_copy_assignment = std::true_type;
        using propagate_on_container_move_assignment = std::true_type;
        using propagate_on_container_swap = std::true_type;
        using is_always_equal = std::true_type;

        template <typename U>
        struct rebind { using other = CuMAllocator<U>; };

        CuMAllocator() noexcept = default;
        template <typename U> constexpr CuMAllocator(const CuMAllocator<U>&) noexcept {}

        T* allocate(std::size_t n) {
            if (n > std::numeric_limits<size_type>::max() / sizeof(T))
                throw std::bad_array_new_length();

            T* p = nullptr;
            cudaError_t err = cudaMallocManaged(&p, n * sizeof(T));
            if (err != cudaSuccess) {
                std::cerr << "cudaMallocManaged failed: " << cudaGetErrorString(err) << '\n';
                throw std::bad_alloc();
            }

            report(p, n, true);
            return p;
        }

        void deallocate(T* p, std::size_t n) noexcept {
            report(p, n, false);
            cudaFree(p);
        }

#if !defined(__CUDA_ARCH__) && __cplusplus < 202002L
        template <typename U, typename... Args>
        void construct(U* p, Args&&... args) {
            ::new((void*)p) U(std::forward<Args>(args)...);
        }

        template <typename U>
        void destroy(U* p) {
            p->~U();
        }
#endif

    private:
        void report(T* p, std::size_t n, bool alloc) const {
#ifndef NDEBUG
            std::cout << (alloc ? "Alloc: " : "Dealloc: ")
                    << n * sizeof(T) << " bytes at " << static_cast<void*>(p) << '\n';
#endif
        }
    };

    template <typename T, typename U>
    bool operator==(const CuMAllocator<T>&, const CuMAllocator<U>&) { return true; }

    template <typename T, typename U>
    bool operator!=(const CuMAllocator<T>&, const CuMAllocator<U>&) { return false; }

} // namespace gpu
