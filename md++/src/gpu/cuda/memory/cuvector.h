/**
 * @file cuvector.h
 * @author poliak
 * CUDA transparent version of std::vector
 */

#pragma once



#include "cumallocator.h"
#include "../../../math/gmath.h"

namespace gpu
{
    template <typename T>
    using cuvector = std::vector< T, gpu::CuMAllocator<T> >;
};

namespace math
{
    /**
     * CUDA allocated managed VArray (double precision)
     */
    typedef VArrayT<Vec, gpu::CuMAllocator > CuDVArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    typedef VArrayT<Vecf, gpu::CuMAllocator > CuFVArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    typedef VArrayT<float3, gpu::CuMAllocator > CuF3VArray;

    /**
     * CUDA allocated managed VArray (single precision)
     */
    typedef VArrayT<double3, gpu::CuMAllocator > CuD3VArray;
}
