/**
 * @file float9.h
 * @author capraz, poliak
 * matrix type single and double precision
 */

#pragma once

#ifndef HOSTDEVICE
  #error "Don't include float9.h without defining HOSTDEVICE"
#else

#define ALIGN16 alignas(16)

#include <type_traits>
/**
 * @struct R9 a matrix of nine floating points
 */
template <typename FP, typename std::enable_if< // allow only floating point types
                                        std::is_floating_point<FP>::value, bool>::type = true>
struct ALIGN16 FP9 {
    /**
     * @brief union allows us to access as xx, xy, xz, ... or (i,j)
     * 
     */
    union {
        struct { FP xx, xy, xz,
                      yx, yy, yz,
                      zx, zy, zz; };
        FP m[3][3];   // row-major 3×3
    };

    /**
     * @brief row getter
     * 
     * @param i row index
     * @return float3 or double3 row copy
     */
    HOSTDEVICE typename std::conditional<
        std::is_same<FP,float>::value, float3, double3>::type
    operator()(int i) const {
        return { m[i][0], m[i][1], m[i][2] };
    }

    // --- element getter (read-only) ---
    /**
     * @brief element getter
     * 
     * @param i row index
     * @param j col index
     * @return float or double element copy
     */
    HOSTDEVICE FP operator()(int i, int j) const {
        return m[i][j];
    }

    /**
     * @brief element reference
     * 
     * @param i row index
     * @param j col index
     * @return float or double element reference
     */
    HOSTDEVICE FP& operator()(int i, int j) {
        return m[i][j];
    }

    /**
     * @brief row setter
     * 
     * @param i row index
     * @param v new row
     */
    HOSTDEVICE void set_row(int i, typename std::conditional<
                                       std::is_same<FP,float>::value, float3, double3>::type const &v) {
        m[i][0] = v.x;
        m[i][1] = v.y;
        m[i][2] = v.z;
    };
};

using float9 = FP9<float>;
using double9 = FP9<double>;

template <typename T9, typename std::enable_if< // allow only float9 and double9
                                        std::is_same<T9,float9>::value ||
                                        std::is_same<T9,double9>::value
                                            , bool>::type = true>
HOSTDEVICE volatile T9& operator+=(volatile T9& a, volatile const T9& b)
{
    a.xx += b.xx;
    a.xy += b.xy;
    a.xz += b.xz;
    a.yx += b.yx;
    a.yy += b.yy;
    a.yz += b.yz;
    a.zx += b.zx;
    a.zy += b.zy;
    a.zz += b.zz;
    return a;
}

// non-volatile variant
template <typename T9, typename std::enable_if< // allow only float9 and double9
                                        std::is_same<T9,float9>::value ||
                                        std::is_same<T9,double9>::value
                                            , bool>::type = true>
HOSTDEVICE T9& operator+=(T9& a, const T9& b)
{
    a.xx += b.xx;
    a.xy += b.xy;
    a.xz += b.xz;
    a.yx += b.yx;
    a.yy += b.yy;
    a.yz += b.yz;
    a.zx += b.zx;
    a.zy += b.zy;
    a.zz += b.zz;
    return a;
}

template <typename T9, typename std::enable_if< // allow only float9 and double9
                                        std::is_same<T9,float9>::value ||
                                        std::is_same<T9,double9>::value
                                            , bool>::type = true>
HOSTDEVICE T9 operator+(T9& a, const T9& b) {
    T9 c;
    c.xx = a.xx + b.xx;
    c.xy = a.xy + b.xy;
    c.xz = a.xz + b.xz;
    c.yx = a.yx + b.yx;
    c.yy = a.yy + b.yy;
    c.yz = a.yz + b.yz;
    c.zx = a.zx + b.zx;
    c.zy = a.zy + b.zy;
    c.zz = a.zz + b.zz;
    return c;
}
#endif