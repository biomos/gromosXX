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
    template <typename T>
    using cuvector = std::vector< T, CuMAllocator<T> >;
    template <typename T>
    using cuhvector = std::vector< T, CuHAllocator<T> >;
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
