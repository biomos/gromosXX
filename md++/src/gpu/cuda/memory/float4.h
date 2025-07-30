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
 * @file float4.h 
 * 4D vector operations
 */

#pragma once

#ifndef HOSTDEVICE
  #error "Don't include float4.h without defining HOSTDEVICE"
#else

/**
 * calculates the scalar (dot) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T dot(const T4 & a, const T4 & b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/**
 * operator to subtract two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return a-b
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator-(const T4 & a, const T4 & b) {
  return T4{
    a.x - b.x,
    a.y - b.y,
    a.z - b.z,
    a.w - b.w
  };
}

/**
 * operator to add two vectors
 * @param[in] a first vector
 * @param[out] b second vector
 * @return a+b
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator+(const T4 & a, const T4 & b) {
  return T4{
    a.x + b.x,
    a.y + b.y,
    a.z + b.z,
    a.w + b.w
  };
}

/**
 * negates a vector
 * @param[in] a the vector
 * @return -a
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator-(const T4 & a) {
  return T4{
    -a.x,
    -a.y,
    -a.z,
    -a.w
  };
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator*(const T4 & a, T b) {
  return T4{
    a.x * b,
    a.y * b,
    a.z * b,
    a.w * b
  };
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator*(T b, const T4 & a) {
  return a * b;
}

/** 
 * scales a vector (division by a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a/b
 */
template <typename T4, typename T = decltype(T4().x),
          typename std::enable_if<
              std::is_same<T4, float4>::value ||
              std::is_same<T4, double4>::value,
              bool>::type = true>
HOSTDEVICE T4 operator/(const T4 & a, T b) {
  b = 1.0f / b;
  return a * b;
}
#endif

