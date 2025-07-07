/**
 * @file cuvector.h
 * @author poliak
 * CUDA transparent version of std::vector
 */

#ifndef INCLUDED_CUVECTOR_H
#define INCLUDED_CUVECTOR_H

#include "cumallocator.h"
#include "math/gmath.h"

namespace gpu
{
    template <typename T>
    using cuvector = std::vector< T, gpu::CuMallocator<T> >;
};

namespace math
{
    /**
     * CUDA allocated VArray (double precision)
     */
    typedef VArrayT<Vec, gpu::CuMallocator > CuDVArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<Vecf, gpu::CuMallocator > CuFVArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<float3, gpu::CuMallocator > CuF3VArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<double3, gpu::CuMallocator > CuD3VArray;
}

#endif
