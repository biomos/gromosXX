/**
 * @file cumallocator.h
 * @author poliak
 * allocator for CUDA managed unified memory, CUDA transparent version of std::vector
 */

#pragma once

namespace cukernel
{  
    template <typename T>
    struct CuMallocator
    {
        typedef T value_type;
    
        CuMallocator () = default;
    
        template<class U>
        constexpr CuMallocator (const CuMallocator <U>&) noexcept {}
    
        T* allocate(std::size_t n)
        {
            if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
                throw std::bad_array_new_length();
    
            T* p = nullptr;
            const cudaError_t error = cudaMallocManaged(&p, n*sizeof(T));
            if (error != cudaSuccess) {
                fprintf(stderr, "code: %d, reason: %s\n", error,
                        cudaGetErrorString(error));
                throw std::bad_alloc();
            }
            report(p, n);
            return p;
        }
    
        void deallocate(T* p, std::size_t n) noexcept
        {
            report(p, n, 0);
            cudaFree(p);
            p = nullptr;
        }
    private:
        void report(T* p, std::size_t n, bool alloc = true) const
        {
        #ifndef NDEBUG
            std::cout << "CuMallocator::" << (alloc ? "Alloc: " : "Dealloc: ") << sizeof(T) * n
                    << " bytes at " << std::hex << std::showbase
                    << reinterpret_cast<void*>(p) << std::dec << '\n';
        #endif
        }
    };

    template<class T, class U>
    bool operator==(const CuMallocator <T>&, const CuMallocator <U>&) { return true; }
    
    template<class T, class U>
    bool operator!=(const CuMallocator <T>&, const CuMallocator <U>&) { return false; }

    template <typename T>
    using cuvector = std::vector< T, cukernel::CuMallocator<T> >;
};
