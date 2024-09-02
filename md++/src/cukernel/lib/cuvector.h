/**
 * @file cuvector.h
 * @author poliak
 * CUDA transparent version of std::vector
 */

#ifndef INCLUDED_CUVECTOR_H
#define INCLUDED_CUVECTOR_H

#include "cumallocator.h"
#include "math/gmath.h"

namespace cukernel
{
    template <typename T>
    using cuvector = std::vector< T, cukernel::CuMallocator<T> >;
};

namespace math
{
    /**
     * CUDA allocated VArray (double precision)
     */
    typedef VArrayT<Vec, cukernel::CuMallocator > CuDVArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<Vecf, cukernel::CuMallocator > CuFVArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<float3, cukernel::CuMallocator > CuFVArray;

    /**
     * CUDA allocated VArray (single precision)
     */
    typedef VArrayT<double3, cukernel::CuMallocator > CuFVArray;
}

#endif
