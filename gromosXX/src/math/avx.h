/**
 * @file avx.h
 * Library for vectorized operations using Intel intrinsics.
 */

#ifndef INCLUDED_AVX_H
#define INCLUDED_AVX_H
#ifdef __AVX2__
#include <immintrin.h>

namespace avx
{
  enum class mode {
    write, add
  };
  const __m256i vec_mask = {-1, -1, -1, 0};

  inline __attribute__((always_inline)) __m256d abs2(__m256d v[3]) {
      //__m256d d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      __m256d d2 =  _mm256_mul_pd(v[0], v[0]);
#ifdef __FMA__
      d2 = _mm256_fmadd_pd(v[1], v[1], d2);
      d2 = _mm256_fmadd_pd(v[2], v[2], d2);
#else
      __m256d tmp = _mm256_mul_pd(v[1], v[1]);
      d2 = _mm256_add_pd(d2, tmp);
      tmp = _mm256_mul_pd(v[2], v[2]);
      d2 = _mm256_add_pd(d2, tmp);
#endif // __FMA__
    return d2;
  }

  /*
    multiple vector orderings are used across the code depending on
    purpose. Since simulation vectors are usually 3D and AVX uses
    2,4 or 8-vectors, to use them efficiently, they sometimes need
    to be packed or unpacked. E.g. for AVX2, 4-vectors are used.
    100% register use is achieved with a set of three 4-vectors,
    with total of 12 doubles. Typical ordering used in the code:
    1. 4 unpacked vectors              - (x0, y0, z0, 0)  ... (x3, y3, z3, 0)
        abbr. U
    2. 3 component-wise packed vectors - (x0, x1, x2, x3) ... (z0, z1, z2, z3)
        abbr. CP
    3. 3 index-wise packed vectors     - (x0, y0, z0, x1) ... (z2, x3, y3, z3)
        abbr. IP
  */

  /**
   * pack 4 U vectors: (x0, y0, z0, 0)  ... (x3, y3, z3, 0)
   * to 3 CP vectors:  (x0, x1, x2, x3) ... (z0, z1, z2, z3)
   * 
   * 4x _mm256_shuffle_pd       LAT 1 THR 1
   * 3x _mm256_permute2f128_pd  LAT 3 THR 1
   * 
   * Example: load positions/forces to process with AVX2
   */
  inline __attribute__((always_inline)) void packCP(const __m256d &a
                                                  , const __m256d &b
                                                  , const __m256d &c
                                                  , const __m256d &d
                                                  , __m256d &dst1
                                                  , __m256d &dst2
                                                  , __m256d &dst3) {
    const __m256d ab1 = _mm256_shuffle_pd(a, b, 0x0);
    const __m256d ab2 = _mm256_shuffle_pd(a, b, 0xF);
    const __m256d cd1 = _mm256_shuffle_pd(c, d, 0x0);
    const __m256d cd2 = _mm256_shuffle_pd(c, d, 0xF);
    dst1 = _mm256_permute2f128_pd(ab1, cd1, 0x20);
    dst2 = _mm256_permute2f128_pd(ab2, cd2, 0x20);
    dst3 = _mm256_permute2f128_pd(ab1, cd1, 0x31);
    // if one needs also the 4th vector
    //dst4 = _mm256_permute2f128_pd(ab2, cd2, 0x31);
  }

  inline __attribute__((always_inline)) void packCP(const __m256d a[4]
                                                  , __m256d dst[3]) {
    packCP(a[0], a[1], a[2], a[3], dst[0], dst[1], dst[2]);
  }

  /**
   * unpack 3 CP vectors: (x0, x1, x2, x3) ... (z0, z1, z2, z3)
   * to 4 U vectors:      (x0, y0, z0, 0)  ... (x3, y3, z3, 0)
   * 
   * 1x _mm256_permute_pd       LAT 1 THR 1
   * 2x _mm256_shuffle_pd       LAT 1 THR 1
   * 4x _mm256_permute2f128_pd  LAT 3 THR 1
   * 
   * Example: store positions/forces to process with AVX2
   */
  inline __attribute__((always_inline)) void unpackCP(const __m256d &x
                                                    , const __m256d &y
                                                    , const __m256d &z
                                                    , __m256d &dst1
                                                    , __m256d &dst2
                                                    , __m256d &dst3
                                                    , __m256d &dst4) {
    const __m256d xy1 = _mm256_shuffle_pd(x, y, 0xC);
    const __m256d xy2 = _mm256_shuffle_pd(x, y, 0x3);
    const __m256d z2 = _mm256_permute_pd(z, 0x5);

    dst1 = _mm256_permute2f128_pd(xy1, z, 0x20);
    dst2 = _mm256_permute2f128_pd(xy2, z2, 0x20);
    dst3 = _mm256_permute2f128_pd(xy2, z, 0x31);
    dst4 = _mm256_permute2f128_pd(xy1, z2, 0x31);
    
    // mask
    const __m256d m256d_mask = _mm256_castsi256_pd(vec_mask);
    dst1 = _mm256_and_pd(dst1, m256d_mask);
    dst2 = _mm256_and_pd(dst2, m256d_mask);
    dst3 = _mm256_and_pd(dst3, m256d_mask);
    dst4 = _mm256_and_pd(dst4, m256d_mask);
  }

  inline __attribute__((always_inline)) void unpackCP(const __m256d v[3]
                                                    , __m256d dst[4]) {
    unpackCP(v[0], v[1], v[2], dst[0], dst[1], dst[2], dst[3]);
  }

  /**
   * unpack 3 IP vectors: (x0, y0, z0, x1) ... (z2, x3, y3, z3)
   * to 4 U vectors:      (x0, y0, z0, 0)  ... (x3, y3, z3, 0)
   * 
   * 1x _mm256_permute_pd       LAT 1 THR 1
   * 2x _mm256_shuffle_pd       LAT 1 THR 1
   * 4x _mm256_permute2f128_pd  LAT 3 THR 1
   * 
   * Example: store positions/forces to process with AVX2
   */
  inline __attribute__((always_inline)) void unpackIP(const __m256d &v1
                                                    , const __m256d &v2
                                                    , const __m256d &v3
                                                    , __m256d &dst1
                                                    , __m256d &dst2
                                                    , __m256d &dst3
                                                    , __m256d &dst4) {
    const __m256d m256d_mask = _mm256_castsi256_pd(vec_mask);

    dst1 = _mm256_and_pd(v1, m256d_mask);

    dst2 = _mm256_blend_pd(v1, v2, 0x7);
    dst2 = _mm256_permute4x64_pd(dst2, 0x93);
    dst2 = _mm256_and_pd(dst2, m256d_mask);

    dst3 = _mm256_permute2f128_pd(v2, v3, 0x21);
    dst3 = _mm256_and_pd(dst3, m256d_mask);

    dst4 = _mm256_permute4x64_pd(v3, 0x39);
    dst4 = _mm256_and_pd(dst4, m256d_mask);
  }

  inline __attribute__((always_inline)) void unpackIP(const __m256d v[3]
                                                    , __m256d dst[4]) {
    unpackIP(v[0], v[1], v[2], dst[0], dst[1], dst[2], dst[3]);
  }



  /**
   * transpose CP : (x0, x1, x2, x3) ... (z0, z1, z2, z3)
   *        to IP : (x0, y0, z0, x1) ... (z2, x3, y3, z3)
   * 
   * 6x _mm256_permute2f128_pd LAT 3 THR 1
   * 7x _mm256_permute4x64_pd LAT 3 THR 1
   * throughput in 13 clocks
   * 
   * Usage: nim for 4 vectors at once
   */
  inline __attribute__((always_inline)) void transposeCP(__m256d &a, __m256d &b, __m256d &c) {
  
    b = _mm256_permute4x64_pd(b, 0xC9);

    __m256d ab = _mm256_permute2f128_pd(a, b, 0x31);
    __m256d ac = _mm256_permute2f128_pd(a, c, 0x21);
    __m256d bc = _mm256_permute2f128_pd(b, c, 0x21);

    ac = _mm256_permute4x64_pd(ac, 0x63);
    ac = _mm256_permute2f128_pd(ac, b, 0x02);
    b = _mm256_permute4x64_pd(ac, 0x78);

    ab = _mm256_permute4x64_pd(ab, 0x2D);
    ab = _mm256_permute2f128_pd(ab, c, 0x30);
    c = _mm256_permute4x64_pd(ab, 0xD2);

    bc = _mm256_permute4x64_pd(bc, 0x78);
    bc = _mm256_permute2f128_pd(bc, a, 0x02);
    a = _mm256_permute4x64_pd(bc, 0x78);
  }
  
  /** transpose CP to IP with copy
   */
  inline __attribute__((always_inline)) void transposeCP_copy(const __m256d v[3], __m256d dst[3]) {
  
    dst[1] = _mm256_permute4x64_pd(v[1], 0xC9);

    __m256d ab = _mm256_permute2f128_pd(v[0], dst[1], 0x31);
    __m256d ac = _mm256_permute2f128_pd(v[0], v[2], 0x21);
    __m256d bc = _mm256_permute2f128_pd(v[1], v[2], 0x21);

    ac = _mm256_permute4x64_pd(ac, 0x63);
    ac = _mm256_permute2f128_pd(ac, dst[1], 0x02);
    dst[1] = _mm256_permute4x64_pd(ac, 0x78);

    ab = _mm256_permute4x64_pd(ab, 0x2D);
    ab = _mm256_permute2f128_pd(ab, v[2], 0x30);
    dst[2] = _mm256_permute4x64_pd(ab, 0xD2);

    bc = _mm256_permute4x64_pd(bc, 0x78);
    bc = _mm256_permute2f128_pd(bc, v[0], 0x02);
    dst[0] = _mm256_permute4x64_pd(bc, 0x78);
  }

  /**
   * outer product and sum from:
   *  IP vector v1: x10, y10, z10, x11, ... z13, x14, y14, z14
   *  IP vector v2: x20, y20, z20, x21, ... z23, x24, y24, z24
   * 
   * to: / x10*x20, x10*y20, x10*z20 \         / x14*x24, x14*y14, x14*z24 \
   *     | y10*x20, y10*y20, y10*z20 | + ... + | y14*x24, y14*y14, y14*z24 |
   *     \ z10*x20, z10*y20, z10*z20 /         \ z14*x24, z14*y14, z14*z24 /
   * in row major order
   * 
   * 20x permute4x64
   * 11x permute2f128
   * 
   * 
   * Usage: virial_tensor
   */
  template <mode M>
  inline __attribute__((always_inline)) void outerIP(const __m256d v1[3], const __m256d v2[3], double mat[9]) {
    static_assert((M == mode::add) || (M == mode::write), "not implemented");
    //static const __m256i mask = {-1,0,0,0};
    // reorder data for multiplication - we get four 3x3 matrices
    __m256d tmp10 = _mm256_permute4x64_pd(v1[0], 0x40);
    __m256d tmp20 = _mm256_permute4x64_pd(v2[0], 0x24);
    //__m256d mat1 = _mm256_mul_pd(tmp10, tmp20);

    __m256d tmp11 = _mm256_permute4x64_pd(v1[0], 0xA5);
    __m256d tmp21 = _mm256_permute4x64_pd(v2[0], 0x49);
    //__m256d mat2 = _mm256_mul_pd(tmp11, tmp21);

    __m256d tmp12 = _mm256_permute4x64_pd(v1[0], 0xFE);
    __m256d tmp22 = _mm256_permute2f128_pd(v2[0], v2[1], 0x21);
    //__m256d mat3 = _mm256_mul_pd(tmp12, tmp22);

    __m256d tmp13 = _mm256_permute4x64_pd(v1[1], 0x40);
    __m256d tmp23 = tmp22;
    tmp23 = _mm256_permute4x64_pd(tmp23, 0x79);
    //__m256d mat4 = _mm256_mul_pd(tmp13, tmp23);

    __m256d tmp14 = _mm256_permute4x64_pd(v1[1], 0xA5);
    __m256d mat5 = _mm256_mul_pd(tmp14, v2[1]);

    __m256d tmp15 = _mm256_permute4x64_pd(v1[1], 0xFE);
    __m256d tmp25 = _mm256_permute2f128_pd(v2[1], v2[2], 0x21);
    __m256d tmp26 = tmp25;
    tmp25 = _mm256_permute4x64_pd(tmp25, 0x92);
    __m256d mat6 = _mm256_mul_pd(tmp15, tmp25);

    __m256d tmp16 = _mm256_permute4x64_pd(v1[2], 0x40);
    __m256d mat7 = _mm256_mul_pd(tmp16, tmp26);

    __m256d tmp17 = _mm256_permute4x64_pd(v1[2], 0xA5);
    __m256d tmp27 = _mm256_permute4x64_pd(v2[2], 0x9E);
    __m256d mat8 = _mm256_mul_pd(tmp17, tmp27);

    __m256d tmp18 = _mm256_permute4x64_pd(v1[2], 0xFE);
    __m256d tmp28 = _mm256_permute4x64_pd(v2[2], 0xE7);
    __m256d mat9 = _mm256_mul_pd(tmp18, tmp28);
    
    // sum up the matrices
    __m256d ires11 = _mm256_permute2f128_pd(mat5, mat6, 0x21);
    __m256d ires12 = _mm256_permute2f128_pd(mat6, mat7, 0x21);
    //ires11 = _mm256_add_pd(ires11, mat1);
    ires11 = _mm256_fmadd_pd(tmp10, tmp20, ires11);
    //ires12 = _mm256_add_pd(ires12, mat2);
    ires12 = _mm256_fmadd_pd(tmp11, tmp21, ires12);


    __m256d ires21 = _mm256_permute2f128_pd(mat7, mat8, 0x21);
    __m256d ires22 = _mm256_permute2f128_pd(mat8, mat9, 0x21);
    //ires21 = _mm256_add_pd(ires21, mat3);
    ires21 = _mm256_fmadd_pd(tmp12, tmp22, ires21);
    //ires22 = _mm256_add_pd(ires22, mat4);
    ires22 = _mm256_fmadd_pd(tmp13, tmp23, ires22);

    __m256d ires3 = _mm256_permute2f128_pd(mat9, ires11, 0x21);
    //ires3 = _mm256_add_pd(ires3, mat5);
    ires3 = _mm256_fmadd_pd(tmp14, v2[1], ires3);

    
    __m256d tmp_res21 = _mm256_permute2f128_pd(ires21, ires22, 0x20);
    tmp_res21 = _mm256_permute4x64_pd(tmp_res21, 0x06);
    tmp_res21 = _mm256_permute2f128_pd(ires21, tmp_res21, 0x21);
    tmp_res21 = _mm256_permute4x64_pd(tmp_res21, 0x93);

    __m256d res1 = _mm256_add_pd(ires11, tmp_res21);

    __m256d tmp_res22 = _mm256_permute2f128_pd(ires22, ires3, 0x20);
    tmp_res22 = _mm256_permute4x64_pd(tmp_res22, 0x06);
    tmp_res22 = _mm256_permute2f128_pd(ires22, tmp_res22, 0x21);
    tmp_res22 = _mm256_permute4x64_pd(tmp_res22, 0x93);

    __m256d res2 = _mm256_add_pd(ires12, tmp_res22);


    __m256d res3 = _mm256_permute4x64_pd(ires3, 0x01);
    res3 = _mm256_add_pd(res3, ires21);

    // add-and-write or just write
    if (M == mode::add) {
      __m256d mat_old1 = _mm256_loadu_pd(&mat[0]);
      __m256d mat_old2 = _mm256_loadu_pd(&mat[4]);
      __m256d mat_old3 = _mm256_broadcast_sd(&mat[8]);

      res1 = _mm256_add_pd(res1, mat_old1);
      res2 = _mm256_add_pd(res2, mat_old2);
      res3 = _mm256_add_pd(res3, mat_old3);
    }

    _mm256_storeu_pd(&mat[0], res1);
    _mm256_storeu_pd(&mat[4], res2);
    _mm_store_sd(&mat[8], _mm256_castpd256_pd128(res3));
  }
} // math
#endif // __AVX2__
#endif // INCLUDED_AVX_H
