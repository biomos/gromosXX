/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file float3.h 
 * 3D vector operations
 */

#pragma once


#ifndef HOSTDEVICE
  #error "Don't include float3.h without defining HOSTDEVICE"
#else

// #include <type_traits>

/**
 * calculates the nearest image for cubic periodic boundary conditions
 * @param[in] i first position
 * @param[in] j second position
 * @param[in] box parameters as float3. Where x is half box, y is the box length, and z is the inverted box length
 * @return the nearest image distance
 */
template <typename T3>
HOSTDEVICE T3 nearest_image(const T3 & i, const T3 & j, const T3 & box_param_x, const T3 & box_param_y, const T3 & box_param_z) {
  T3 r;

  r.x = i.x - j.x;
  if (fabs(r.x) > box_param_x.x)
    r.x -= box_param_x.y * rintf(r.x * box_param_x.z);

  r.y = i.y - j.y;
  if (fabs(r.y) > box_param_y.x)
    r.y -= box_param_y.y * rintf(r.y * box_param_y.z);

  r.z = i.z - j.z;
  if (fabs(r.z) > box_param_z.x)
    r.z -= box_param_z.y * rintf(r.z * box_param_z.z);

  return r;
}

// HOSTDEVICE float3 nearest_image(const float3 &a, const float3 &b) {
//   float3 nim;
//   nim.x = a.x - b.x;
//   nim.y = a.y - b.y;
//   nim.z = a.z - b.z;
//   float* fnim = reinterpret_cast<float*>(&nim);
//   float* box_half = reinterpret_cast<float*>(&(device_param.box.half));
//   float* box_full = reinterpret_cast<float*>(&(device_param.box.full));
//   #ifdef NDEBUG
//   #pragma unroll(3)
//   #endif
//   for (unsigned i = 0; i < 3; ++i) {
//     while (fnim[i] > box_half[i]) {
//       fnim[i] -= box_full[i];
//     }
//     while (fnim[i] < -box_half[i]) {
//       fnim[i] += box_full[i];
//     }
//   }
//   return nim;
// }

/**
 * calculates the scalar (dot) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE auto dot(const T3 & a, const T3 & b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * calculates the vector (cross) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the cross product
 */
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T3 cross(const T3 & a, const T3 & b) {
  return T3{
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  };
}

/**
 * operator to subtract two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return a-b
 */
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator-(const T3 & a, const T3 & b) {
  return T3{
    a.x - b.x,
    a.y - b.y,
    a.z - b.z
  };
}

/**
 * operator to add two vectors
 * @param[in] a first vector
 * @param[out] b second vector
 * @return a+b
 */
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator+(const T3 & a, const T3 & b) {
  return T3{
    a.x + b.x,
    a.y + b.y,
    a.z + b.z
  };
}

/**
 * negates a vector
 * @param[in] a the vector
 * @return -a
 */
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator-(const T3 & a) {
  return T3{
    -a.x,
    -a.y,
    -a.z
  };
}

// non-volatile variant
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                                (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
                                                (std::is_same<T,float>::value  || std::is_same<T,double>::value)
                                            , bool>::type = true>
HOSTDEVICE T3& operator-(const T3& a, T b) {
  return T3{
    a.x-b,
    a.y-b,
    a.z-b
  };
}

// non-volatile variant
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                                (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
                                                (std::is_same<T,float>::value  || std::is_same<T,double>::value)
                                            , bool>::type = true>
HOSTDEVICE T3& operator-(T b, const T3& a) {
  return T3{
    b-a.x,
    b-a.y,
    b-a.z
  };
}

// non-volatile variant
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                                (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
                                                (std::is_same<T,float>::value  || std::is_same<T,double>::value)
                                            , bool>::type = true>
HOSTDEVICE T3& operator-=(T3& a, T b) {
    a.x -= b;
    a.y -= b;
    a.z -= b;
    return a;
}

// non-volatile variant
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                  std::is_same<T3,float3>::value || std::is_same<T3,double3>::value
                              , bool>::type = true>
HOSTDEVICE T3& operator-=(T3& a, const T3& b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
template <typename T3, typename T,
          typename CompT = decltype(T3{}.x),
          typename std::enable_if<
              (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
              std::is_same<T, CompT>::value,
              bool>::type = true>
HOSTDEVICE T3 operator*(const T3 & a, T b) {
  return T3{
    a.x * b,
    a.y * b,
    a.z * b
  };
}

// non-volatile variant
template <typename T3, typename T,
          typename CompT = decltype(T3{}.x),
          typename std::enable_if<
              (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
              std::is_same<T, CompT>::value,
              bool>::type = true>
HOSTDEVICE T3& operator*=(T3& a, T b) {
    a.x *= b;
    a.y *= b;
    a.z *= b;
    return a;
}

template <typename T3, typename T,
          typename CompT = decltype(T3{}.x),
          typename std::enable_if<
              (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
              std::is_same<T, CompT>::value,
              bool>::type = true>
HOSTDEVICE T3 operator*(T b, const T3 & a) {
    return a * b;
}

/**
 * addition operators for float3 and double3 - volatile variant, uncomment only if needed
 */
// template <typename T3, typename std::enable_if< // allow only float3 and double3
//                                         std::is_same<T3,float3>::value ||
//                                         std::is_same<T3,double3>::value
//                                             , bool>::type = true>
// HOSTDEVICE volatile T3& operator+=(volatile T3& a, volatile const T3& b) {
//     a.x += b.x;
//     a.y += b.y;
//     a.z += b.z;
//     return a;
// }

// non-volatile variant
template <typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T3& operator+=(T3& a, const T3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

// // volatile variant
// template <typename T, typename T3, typename std::enable_if< ... , bool>::type = true>
// HOSTDEVICE volatile T3& operator*=(volatile T3& a, T b) {
//     a.x *= b;
//     a.y *= b;
//     a.z *= b;
//     return a;
// }

// Helper: compute reciprocal of b in the same type as T3's component
template <typename T, typename T3>
HOSTDEVICE auto reciprocal(T b) -> decltype(T3{}.x) {
    using CompT = decltype(T3{}.x); // float for float3, double for double3
    CompT rcp;

#ifdef __CUDA_ARCH__
    if constexpr (std::is_same<CompT,float>::value) {
        rcp = __frcp_rn(static_cast<float>(b)); // fast float reciprocal on device
    } else {
        rcp = static_cast<CompT>(1) / static_cast<CompT>(b);
    }
#else
    rcp = static_cast<CompT>(1) / static_cast<CompT>(b); // host division
#endif

    return rcp;
}

/** 
 * scales a vector (division by a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a/b
 */
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                                (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
                                                std::is_arithmetic<T>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator/(const T3& a, T b) {
    using CompT = decltype(T3{}.x);
    // CompT rcp = static_cast<CompT>(1.) / b;
    CompT rcp = reciprocal<T, T3>(b);
    return a * rcp;
}

// non-volatile variant
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                                (std::is_same<T3,float3>::value || std::is_same<T3,double3>::value) &&
                                                std::is_arithmetic<T>::value
                                            , bool>::type = true>
HOSTDEVICE T3& operator/=(T3& a, T b) {
    using CompT = decltype(T3{}.x);
    // CompT rcp = static_cast<CompT>(1.) / b;
    CompT rcp = reciprocal<T, T3>(b);
    a.x *= rcp;
    a.y *= rcp;
    a.z *= rcp;
    return a;
}

// // Helper traits to simplify casting to T3
// template<typename, typename = void>
// struct has_xyz_only : std::false_type {};

// template<typename T>
// struct has_xyz_only<T, std::void_t<
//     decltype(std::declval<T>().x),
//     decltype(std::declval<T>().y),
//     decltype(std::declval<T>().z)
// >> {
//     static constexpr bool has_w = requires(T t) { t.w; };
//     static constexpr bool value = !has_w;
// };

#endif

