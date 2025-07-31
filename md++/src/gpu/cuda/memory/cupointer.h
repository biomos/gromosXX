/**
 * @file cupointer.h
 * @author poliak
 * implementation of smart pointers for CUDA
 * It can do just standard deallocation, or even tracked one, where it calls a memory manager
 * to deallocate tracked memory gracefully
 */

#pragma once

#include <memory>

namespace gpu {
    template <typename T, typename MM = void>
    struct CuDeleter;
    // Deallocate gracefully
    // Case 1: No memory manager
    template <typename T>
    struct CuDeleter<T, void> {
        void operator()(T* ptr) const {
            if (ptr) {
                cudaError_t err = cudaFree(ptr);
                if (err != cudaSuccess) {
                    std::cerr << "cudaFree failed [" << __FILE__ << ":" << __LINE__ << "]" << cudaGetErrorString(err) << std::endl;
                }
            }
        }
    };

    // Case 2: With memory manager
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam MM - Memory manager
     */
    template <typename T, typename MM>
    struct CuDeleter {
        MM* mgr;
        explicit CuDeleter(MM* m = nullptr) : mgr(m) {}

        void operator()(T* ptr) const {
            if (ptr) {
                if (mgr) {
                    mgr->deallocate(ptr);
                } else {
                    cudaError_t err = cudaFree(ptr);
                    if (err != cudaSuccess) {
                        std::cerr << "cudaFree failed [" << __FILE__ << ":" << __LINE__ << "]" << cudaGetErrorString(err) << std::endl;
                    }
                }
            }
        }
    };
    // For single objects
    template <typename T>
    using cu_uniquevar_ptr = std::unique_ptr<T, CuDeleter<T>>;

    // For arrays
    template <typename T>
    using cu_uniquearr_ptr = std::unique_ptr<T[], CuDeleter<T>>;

    template <typename T>
    using cu_shared_ptr = std::shared_ptr<T[]>;

    template <typename T>
    using cu_weak_ptr = std::weak_ptr<T[]>;
}