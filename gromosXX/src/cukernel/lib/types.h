/**
 * @file types.h
 * additional types of general use
 */

#ifndef INCLUDED_TYPES_H
#define INCLUDED_TYPES_H

#include <type_traits>

#define HOSTDEVICE __device__ __forceinline__
#include "float3.h"
#include "float9.h"
#include "double3.h"

template <typename T2, typename std::enable_if< // allow only float2 and double2
                                        std::is_same<T2,float2>::value ||
                                        std::is_same<T2,double2>::value
                                            , bool>::type = true>
HOSTDEVICE T2 operator+(T2& a, const T2& b) {
    T2 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

/**
 * addition operators for float2 and double2
 */
template <typename T2, typename std::enable_if< // allow only float2 and double2
                                        std::is_same<T2,float2>::value ||
                                        std::is_same<T2,double2>::value
                                            , bool>::type = true>
HOSTDEVICE volatile T2& operator+=(volatile T2& a, volatile const T2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

// non-volatile variant
template <typename T2, typename std::enable_if< // allow only float2 and double2
                                        std::is_same<T2,float2>::value ||
                                        std::is_same<T2,double2>::value
                                            , bool>::type = true>
HOSTDEVICE T2& operator+=(T2& a, const T2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}


#undef HOSTDEVICE

#endif /* TYPES_H */
