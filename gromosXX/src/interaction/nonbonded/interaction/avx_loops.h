/**
 * @file avx_loops.h
 */

#ifndef INCLUDED_AVX_LOOPS_H
#define INCLUDED_AVX_LOOPS_H

#ifdef __AVX2__
#include <immintrin.h>
#include "../../../math/avx.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

namespace interaction {
  
  void lj_crf_fast_solute_loop_avx2(topology::Topology &topo,
                                    configuration::Configuration &conf,
                                    simulation::Simulation &sim,
                                    const Pairlist &pairlist,
                                    Storage &storage,
                                    math::Periodicity<math::rectangular> &periodicity,
                                    Nonbonded_Innerloop
                                      <interaction::Interaction_Spec
                                          <math::rectangular, simulation::lj_crf_func>
                                      > &innerloop) {

  std::vector<unsigned>::const_iterator j_it, j_to;

  unsigned size_i = unsigned(pairlist.size());
  DEBUG(10, "lj_crf2 outerloop pairlist size " << size_i);

  const unsigned end = topo.num_solute_atoms();

  const double* pos0 = &conf.current().pos(0)[0];
  constexpr size_t pos_size = sizeof(conf.current().pos(0));
  
  //__m256d m256d_groupForce[27];
  alignas(32) double groupForce[84];

  for (unsigned ai = 0; ai < end; ++ai) {
    // zero groupforce
    for (unsigned i = 0; i < 84; i += 4) {
      _mm256_storeu_pd(groupForce + i, _mm256_setzero_pd());
    }
    const math::Vec &posI = conf.current().pos(ai);
    const unsigned eg_i = topo.atom_energy_group(ai);
    // load position of first atom
    __m256d m256d_ipos[3];
    for (unsigned i = 0; i < 3; ++i) {
      m256d_ipos[i] = _mm256_set1_pd(posI[i]);
    }

    for (j_it = pairlist[ai].begin(), j_to = pairlist[ai].end();
         j_it < j_to; j_it += 4) {
      int remaining = j_to - j_it;

      __m128i m128i_mask;
      m128i_mask = _mm_cmpeq_epi16(m128i_mask, m128i_mask); // fills vector with ones
      __m128i m128i_js;
      __m256d m256d_jpos[3];
      // load mask to 128-bit variable
      if (remaining >= 4) {
        DEBUG(7, "Full block, remaining: " << remaining);
        // load integers of 4 atoms
        m128i_js = _mm_lddqu_si128((__m128i*)&*j_it);
        DEBUG(10, "Loading j atoms: " << ((int*)&m128i_js)[0] << " " << ((int*)&m128i_js)[1] << " " << ((int*)&m128i_js)[2] << " " << ((int*)&m128i_js)[3]);
        // load positions based on integers
        for (unsigned i = 0; i < 3; ++i) {
          m256d_jpos[i] = _mm256_i32gather_pd(pos0 + i, m128i_js * pos_size, 1);
        }
      }
      else {        
        int cnt = 4 - remaining;
        while(cnt--) {
          m128i_mask = _mm_srli_si128(m128i_mask, 4);
        }
        // load integers of the remaining atoms
        const __m128i vindex = _mm_set_epi32(12, 8, 4, 0);
        m128i_js = _mm_mask_i32gather_epi32(_mm_setzero_si128(), (int*)&*j_it, vindex, m128i_mask, 1);
        
        const __m256i m256i_mask = _mm256_cvtepi32_epi64(m128i_mask);
        // load positions based on integers
        for (unsigned i = 0; i < 3; ++i) {
          m256d_jpos[i] = _mm256_mask_i32gather_pd(_mm256_setzero_pd(), pos0 + i, m128i_js * pos_size, _mm256_castsi256_pd(m256i_mask), 1);
        }
      }
      // nearest image vectors 0:(x0 ... x3), 1:(y0 ... y3), 2:(z0 ... z3)
      __m256d m256d_nim[3];

      // shift indices for (atom_j0 ... atom_j3)
      __m128i m128i_k;
      periodicity.nearest_image(m256d_ipos, m256d_jpos, m256d_nim, &m128i_k);

      // distances between atom_i and (atom_j0 ... atom_j3)
      const __m256d m256d_nim_d2 = avx::abs2(m256d_nim);
      DEBUG(15, "m256d_nim_d2: " << m256d_nim_d2[0] << " " << m256d_nim_d2[1] << " " << m256d_nim_d2[2] << " " << m256d_nim_d2[3]);
      
      // force and energies (scalars) between atom_i and (atom_j0 ... atom_j3)
      __m256d m256d_f, m256d_e_lj, m256d_e_crf;

      innerloop.lj_crf_innerloop_avx(topo, ai, m128i_js, m256d_nim_d2, m128i_mask, m256d_f, m256d_e_lj, m256d_e_crf);

      int* const js = (int*)&m128i_js;
      for (unsigned i = 0; i < 4; ++i) {
        DEBUG(15, "interaction " << ai << "-" << js[i]);
        DEBUG(15, "d2:    " << m256d_nim_d2[i]);
        DEBUG(15, "f:     " << m256d_f[i]);
        DEBUG(15, "e_lj:  " << m256d_e_lj[i]);
        DEBUG(15, "e_crf: " << m256d_e_crf[i]);
      }
      // gather atoms energy groups and store energies
      constexpr size_t eg_size = sizeof(topo.atom_energy_group()[0]);
      __m128i m128i_eg_j = _mm_i32gather_epi32((int*)&topo.atom_energy_group()[0], m128i_js * eg_size, 1);
      int* const eg_j = (int*)&m128i_eg_j;
      double* const e_lj = (double*)&m256d_e_lj;
      double* const e_crf = (double*)&m256d_e_crf;
      const int imask = _mm_movemask_epi8(m128i_mask);
      for (unsigned i = 0; i < 4; ++i) {
        if (imask>>i*4&1) {
          storage.energies.lj_energy[eg_i][eg_j[i]] += e_lj[i];
          storage.energies.crf_energy[eg_i][eg_j[i]] += e_crf[i];
        }
      }
      
      // packed force 0:(x0 ... x3), 1:(y0 ... y3), 2:(z0 ... z3)
      __m256d m256d_forceCP[3];
      for (unsigned i = 0; i < 3; ++i) {
        m256d_forceCP[i] = _mm256_mul_pd(m256d_f, m256d_nim[i]);
      }
      int* const js = (int*)&m128i_js;
      for (unsigned i = 0; i < 4; ++i) {
        DEBUG(15, "interaction " << ai << "-" << js[i]);
        DEBUG(15, "nim: " << m256d_nim[0][i] << " " << m256d_nim[1][i] << " "  << m256d_nim[2][i]);
        DEBUG(15, "f: " << m256d_forceCP[0][i] << " " << m256d_forceCP[1][i] << " "  << m256d_forceCP[2][i]);
      }

      // unpacked force i:(xi, yi, zi, junk)
      __m256d m256d_force[4];
      avx::unpackCP(m256d_forceCP, m256d_force);
      m128i_k = _mm_add_epi32(m128i_k, _mm_set1_epi32(13));
      int* const k = (int*)&m128i_k;

      for (unsigned i = 0; i < 4; ++i) {
        if (imask>>i*4&1) {
          // write i-th groupForce
          int* const ki = k + i;
          double* gf = groupForce + *ki * 3;
          __m256d m256d_gf = _mm256_loadu_pd(gf);
          m256d_gf = _mm256_add_pd(m256d_gf, m256d_force[i]);
          _mm256_storeu_pd(gf, m256d_gf);
          DEBUG(15, "avx::vec_mask: " << avx::vec_mask[0] << " " << avx::vec_mask[1] << " " << avx::vec_mask[2] << " " << avx::vec_mask[3]);
          // load stored j force
          int* const jsi = js + i;
          if (*ki != 13) {
            DEBUG(10, "*ki != 13, " << "ai: " << ai << ", aj: " << *jsi);
          }
          double* const f_ptr = &storage.force(*jsi)[0];
          __m256d f = _mm256_loadu_pd(f_ptr);
          // load unaligned, masked add/sub (using AND op) and write back
          // the fourth (pseudo-masked) element is written back unchanged
          f = _mm256_sub_pd(f, m256d_force[i]);
          _mm256_storeu_pd(f_ptr, f);
        }
      }
    }
    for (unsigned i = 0; i < 84; i += 4) {
      DEBUG(15, i << "-" << i + 3 << ":\t" << groupForce[i] << " " << groupForce[i+1] << " " << groupForce[i+2] << " " << groupForce[i+3]);
    }
    // sum up groupForce, write to virial tensor and write to ai-th force
    // reduce over triples of packed transposed forces (12 packed elements)
    // IP packed force 0:(x0, y0, z0, x1), 1:(y1 ... y3), 2:(z3, x4, y4, z4)
    __m256d m256d_forceIP[3];
    // initialize
    for (unsigned i = 0; i < 3; ++i) {
      m256d_forceIP[i] = _mm256_setzero_pd();
    }

    for (unsigned i = 0; i < 7; ++i) {
      __m256d m256d_f[3];
      const double* shift0 = &periodicity.shift(i*4).pos[0];
      constexpr size_t shift_size = sizeof(periodicity.shift(0));
      __m256d m256d_shift[3];
      const __m128i vindex = _mm_set_epi32(3 * shift_size, 2 * shift_size, shift_size, 0);
      for (unsigned j = 0; j < 3; ++j) {
        double* const gf = groupForce + i*12 + j*4;
        m256d_f[j] = _mm256_loadu_pd(gf);
        m256d_forceIP[j] = _mm256_add_pd(m256d_forceIP[j], m256d_f[j]);
        m256d_shift[j] = _mm256_i32gather_pd(shift0 + j, vindex, 1);
      }
      // transpose from x1, x2, x3, x4 ... z1, z2, z3, z4
      //             to x1, y1, z1, x2 ... z3, x4, y4, z4
      avx::transposeCP(m256d_shift[0], m256d_shift[1], m256d_shift[2]);
      math::Matrix* virial_tensor = &storage.virial_tensor;
      avx::outerIP<avx::mode::add>(m256d_shift, m256d_f, (double*)virial_tensor);
    }

    //unpack to individual atomic forces
    __m256d m256d_force[4];
    avx::unpackIP(m256d_forceIP, m256d_force);
    double* const f_ptr = &storage.force(ai)[0];
    __m256d m256d_force_ai = _mm256_loadu_pd(f_ptr);

    // sum them up
    for (unsigned i = 0; i < 4; ++i) {
      // again, fourth value is zeroed, so the fourth element is written back unchanged
      m256d_force_ai = _mm256_add_pd(m256d_force_ai, m256d_force[i]);
    }
    // and write
    _mm256_storeu_pd(f_ptr, m256d_force_ai);
  }
  // now do the standard virial tensor - the shift part was already done
  // we will do 4 atoms at once - 12 coordinates
  // loading 4 doubles at once, 3 vectors of 4 doubles, total 12 doubles
  math::Matrix* virial_tensor = &storage.virial_tensor;
  const unsigned num_iterations = topo.num_solute_atoms() / 4;
  const unsigned last_size = topo.num_solute_atoms() % 4;
  for (unsigned i = 0; i < num_iterations; ++i) {
    __m256d m256d_pos[3];
    __m256d m256d_force[3];
    constexpr size_t vec_size = sizeof(math::Vec);
    static_assert(vec_size == 3*sizeof(double), "the avx loop relies on math::Vec without padding");
    const math::Vec* pos = &conf.current().pos(4*i);
    const math::Vec* force = &storage.force(4*i);
    for (unsigned j = 0; j < 3; ++j) {
      const double* dp = (double*)pos + 4*j;
      const double* df = (double*)force + 4*j;
      m256d_pos[j] = _mm256_loadu_pd(dp);
      m256d_force[j] = _mm256_loadu_pd((double*)df);
    }
    avx::outerIP<avx::mode::add>(m256d_pos, m256d_force, (double*)virial_tensor);
  }
  // last iteration, load and zero non-relevant
  // gather supports unaligned mask load
  constexpr size_t dsize = sizeof(double);
  // initialize vindex
  __m128i m128i_vindex = _mm_set_epi32(3,2,1,0) * dsize;
  __m256d m256d_pos[3];
  __m256d m256d_force[3];
  __m256i m256i_mask[3];
  for (unsigned i = 0; i < 3; ++i)
    m256i_mask[i] = _mm256_setzero_si256();
  long long* imask = (long long*)&m256i_mask;
  for (unsigned i = 0; i < last_size; ++i) {
    imask[i] = -1;
  }
  const math::Vec* pos = &conf.current().pos(num_iterations*4);
  const math::Vec* force = &storage.force(num_iterations*4);
  for(unsigned i = 0; i < 3; ++i) {
    const double* dp = (double*)pos + i*4;
    const double* df = (double*)force + i*4;
    m256d_pos[i] =   _mm256_mask_i32gather_pd(_mm256_setzero_pd(), dp, m128i_vindex, _mm256_castsi256_pd(m256i_mask[i]), 1);
    m256d_force[i] = _mm256_mask_i32gather_pd(_mm256_setzero_pd(), df, m128i_vindex, _mm256_castsi256_pd(m256i_mask[i]), 1);
  }
  avx::outerIP<avx::mode::add>(m256d_pos, m256d_force, (double*)virial_tensor);
  return;
  }
} // interaction
#endif // __AVX2__
#endif // INCLUDED_AVX_H
