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
 * @file morton.h
 * implementation of Morton ordering for efficient memory locality
 */

#include "gpu/cuda/cuhostdevice.h"

/**
 * 32-bit version
 */
// "split by 3" helper: spreads bits of n so there are 2 zeros between each bit
HOSTDEVICE uint32_t split_by_3_32(uint32_t n) {
    /**
     * possibly for CPU just
     * 
    #include <immintrin.h>
    return _pdep_u32(x, 0x09249249u);
     */
    n &= 0x3FFu;                  // keep only 10 bits
    n = (n | (n << 16)) & 0x30000FFu;
    n = (n | (n << 8))  & 0x300F00Fu;
    n = (n | (n << 4))  & 0x30C30C3u;
    n = (n | (n << 2))  & 0x9249249u;
    return n;
}

// Reverse of split_by_3: compact bits by removing 2 zeros between them.
HOSTDEVICE uint32_t compact_by_3_32(uint32_t n) {
    /**
     * possibly for CPU just
     * 
    #include <immintrin.h>
    return _pext_u32(x, 0x09249249u);
     */
    n &= 0x9249249u;
    n = (n ^ (n >> 2))  & 0x30C30C3u;
    n = (n ^ (n >> 4))  & 0x300F00Fu;
    n = (n ^ (n >> 8))  & 0x30000FFu;
    n = (n ^ (n >> 16)) & 0x3FFu;
    return n;
}

// Morton index for 3D cell (i,j,k)
HOSTDEVICE uint32_t morton_index_32(uint32_t i, uint32_t j, uint32_t k) {
    return (split_by_3_32(k) << 2) | (split_by_3_32(j) << 1) | split_by_3_32(i);
}

// Morton decode: return (x,y,z)
HOSTDEVICE void morton_decode_32(uint32_t index, uint32_t &i, uint32_t &j, uint32_t &k) {
    i = compact_by_3_32(index);
    j = compact_by_3_32(index >> 1);
    k = compact_by_3_32(index >> 2);
}

/**
 * 16-bit version
 */
// "split by 3 (16-bit)" helper: spreads bits of n so there are 2 zeros between each bit
HOSTDEVICE uint16_t split_by_3_16(uint16_t n) {
    n &= 0x1Fu;                   // keep only 5 bits (max 32 per dimension)
    n = (n | (n << 8))  & 0x010F; // spread bits
    n = (n | (n << 4))  & 0x0303;
    n = (n | (n << 2))  & 0x09249; 
    return n;
}

// Reverse of split_by_3: compact bits by removing 2 zeros between them.
HOSTDEVICE uint16_t compact_by_3_16(uint16_t n) {
    n &= 0x09249u;
    n = (n ^ (n >> 2))  & 0x0303u;
    n = (n ^ (n >> 4))  & 0x010Fu;
    n = (n ^ (n >> 8))  & 0x1Fu;
    return n;
}

// Morton index for 3D cell (i,j,k) → 16-bit version
HOSTDEVICE uint16_t morton_index_16(uint16_t i, uint16_t j, uint16_t k) {
    return (split_by_3_16(k) << 2) | (split_by_3_16(j) << 1) | split_by_3_16(i);
}

// Morton decode: return (i,j,k) from 16-bit Morton index
HOSTDEVICE void morton_decode_16(uint16_t code, uint16_t &i, uint16_t &j, uint16_t &k) {
    i = compact_by_3_16(code);
    j = compact_by_3_16(code >> 1);
    k = compact_by_3_16(code >> 2);
}