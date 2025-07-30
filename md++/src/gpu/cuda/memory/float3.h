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
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value
                                            , bool>::type = true>
HOSTDEVICE T dot(const T3 & a, const T3 & b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * calculates the vector (cross) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
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

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value &&
                                        std::is_same<T,float>::value ||
                                        std::is_same<T,double>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator*(const T3 & a, T b) {
  return T3{
    a.x * b,
    a.y * b,
    a.z * b
  };
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value &&
                                        std::is_same<T,float>::value ||
                                        std::is_same<T,double>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator*(T b, const T3 & a) {
  return a * b;
}

/** 
 * scales a vector (division by a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a/b
 */
template <typename T, typename T3, typename std::enable_if< // allow only float3 and double3
                                        std::is_same<T3,float3>::value ||
                                        std::is_same<T3,double3>::value &&
                                        std::is_same<T,float>::value ||
                                        std::is_same<T,double>::value
                                            , bool>::type = true>
HOSTDEVICE T3 operator/(const T3 & a, T b) {
  b = 1.0f / b;
  return a * b;
}
#endif

