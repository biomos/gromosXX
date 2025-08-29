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
 * @file box.h
 * A light-weight Box struct for GPU
 * thin version of math::Box
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/types.h"

namespace gpu {
    struct Box {
        using value_type  = FPL_TYPE; // float
        using vec_type    = FP3<FPL_TYPE>; // float3
        using matrix_type = FP9<FPL_TYPE>; // float9

    private:
        matrix_type d_b;

    public:
        // constructors
        HOSTDEVICE Box() : d_b{0,0,0,0,0,0,0,0,0} {}

        HOSTDEVICE explicit Box(value_type d) {
            d_b.xx = d; d_b.xy = d; d_b.xz = d;
            d_b.yx = d; d_b.yy = d; d_b.yz = d;
            d_b.zx = d; d_b.zy = d; d_b.zz = d;
        }

        HOSTDEVICE Box(vec_type const &v1, vec_type const &v2, vec_type const &v3) {
            d_b.set_row(0, v1);
            d_b.set_row(1, v2);
            d_b.set_row(2, v3);
        }

        HOSTDEVICE Box(const Box &other) : d_b(other.d_b) {}

        __host__ inline Box(const math::Box &other) {
            d_b = other;
        }

        // row access
        HOSTDEVICE vec_type operator()(int i) const { return d_b(i); }
        HOSTDEVICE void set_row(int i, const vec_type &v) { d_b.set_row(i, v); }

        // element access
        HOSTDEVICE value_type operator()(int i, int j) const { return d_b(i,j); }
        HOSTDEVICE value_type& operator()(int i, int j) { return d_b(i,j); }

        // arithmetic operators
        HOSTDEVICE Box& operator+=(const Box &other) {
            for(int i=0;i<3;i++){
                vec_type v = d_b(i);
                vec_type u = other.d_b(i);
                d_b.set_row(i, vec_type{v.x+u.x, v.y+u.y, v.z+u.z});
            }
            return *this;
        }

        HOSTDEVICE Box& operator-=(const Box &other) {
            for(int i=0;i<3;i++){
                vec_type v = d_b(i);
                vec_type u = other.d_b(i);
                d_b.set_row(i, vec_type{v.x-u.x, v.y-u.y, v.z-u.z});
            }
            return *this;
        }

        HOSTDEVICE Box& operator*=(value_type s) {
            for(int i=0;i<3;i++){
                vec_type v = d_b(i);
                d_b.set_row(i, vec_type{v.x*s, v.y*s, v.z*s});
            }
            return *this;
        }

        HOSTDEVICE Box& operator/=(value_type s) {
            value_type inv = 1.0/s;
            return (*this *= inv);
        }

        // non-member friend operators
        friend HOSTDEVICE Box operator+(Box lhs, const Box &rhs) { return lhs += rhs; }
        friend HOSTDEVICE Box operator-(Box lhs, const Box &rhs) { return lhs -= rhs; }
        friend HOSTDEVICE Box operator*(Box lhs, value_type s) { return lhs *= s; }
        friend HOSTDEVICE Box operator*(value_type s, Box rhs) { return rhs *= s; }
        friend HOSTDEVICE Box operator/(Box lhs, value_type s) { return lhs /= s; }
    };
} // namespace gpu