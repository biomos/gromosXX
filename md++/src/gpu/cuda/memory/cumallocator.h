/**
 * @file cumallocator.h
 * @author poliak
 * allocator for CUDA managed unified memory, CUDA transparent version of std::vector
 */

#ifndef INCLUDED_CUMALLOCATOR_H
#define INCLUDED_CUMALLOCATOR_H

#include <cuda_runtime.h>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

namespace cuda {

template <typename T>
class CuMallocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::true_type;

    template <typename U>
    struct rebind { using other = CuMallocator<U>; };

    CuMallocator() noexcept = default;
    template <typename U> constexpr CuMallocator(const CuMallocator<U>&) noexcept {}

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

#ifndef __CUDA_ARCH__
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
bool operator==(const CuMallocator<T>&, const CuMallocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const CuMallocator<T>&, const CuMallocator<U>&) { return false; }

template <typename T>
using cuvector = std::vector<T, CuMallocator<T>>;

} // namespace cuda

#endif
