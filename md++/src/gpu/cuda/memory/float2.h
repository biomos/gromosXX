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
 * @file float2.h 
 * 2D vector operations
 */

#pragma once

#ifndef HOSTDEVICE
  #error "Don't include float2.h without defining HOSTDEVICE"
#else

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

#endif

