/**
 * @file float9.h
 * @author capraz, poliak
 * matrix type single and double precision
 */

#ifndef INCLUDED_FLOAT9_H
#define INCLUDED_FLOAT9_H

#ifndef HOSTDEVICE
#error "Don't include float9.h without defining HOSTDEVICE"
#else

#include <type_traits>
/**
 * @struct float9 a matrix in single precision
 */
template <typename R, typename std::enable_if< // allow only floating point types
                                        std::is_floating_point<R>::value, bool>::type = true>
struct R9 {
    float xx;
    float xy;
    float xz;
    float yx;
    float yy;
    float yz;
    float zx;
    float zy;
    float zz;
};

typedef R9<float> float9;
typedef R9<double> double9;

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
#endif