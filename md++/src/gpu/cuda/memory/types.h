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
 * @file types.h
 * additional types of general use
 */

#pragma once

#include <type_traits>

#define HOSTDEVICE __host__ __device__ __forceinline__
#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float9.h"

/**
 * calculates the squared length of a vector
 * @param a the vector
 * @return squared length
 */
template <typename V, typename T = decltype(V().x),
          typename std::enable_if<
              std::is_same<V, float2>::value ||
              std::is_same<V, float3>::value ||
              std::is_same<V, float4>::value ||
              std::is_same<V, double2>::value ||
              std::is_same<V, double3>::value ||
              std::is_same<V, double4>::value,
              bool>::type = true>
HOSTDEVICE T abs2(const V & a) {
  return dot(a,a);
}

/**
 * calculates the length of a vector
 * @param a the vector
 * @return length
 */
template <typename V, typename T = decltype(V().x),
          typename std::enable_if<
              std::is_same<V, float2>::value ||
              std::is_same<V, float3>::value ||
              std::is_same<V, float4>::value ||
              std::is_same<V, double2>::value ||
              std::is_same<V, double3>::value ||
              std::is_same<V, double4>::value,
              bool>::type = true>
HOSTDEVICE T abs(const V & a) {
  return sqrtf(abs2(a));
}

#undef HOSTDEVICE

