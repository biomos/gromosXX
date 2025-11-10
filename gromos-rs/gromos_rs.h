#ifndef GROMOS_RS_H
#define GROMOS_RS_H

#pragma once

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * C FFI: LJ + CRF inner loop
 *
 * # Safety
 * Caller must ensure all pointers are valid and have correct length
 *
 * # Parameters
 * * `positions` - Flat array of positions [x0,y0,z0,x1,y1,z1,...] of length n_atoms*3
 * * `charges` - Array of charges, length n_atoms
 * * `iac` - Integer atom codes (atom types), length n_atoms
 * * `n_atoms` - Number of atoms
 * * `pairlist` - Flat array of pairs [i0,j0,i1,j1,...] of length n_pairs*2
 * * `n_pairs` - Number of pairs
 * * `lj_params` - Flat array of LJ parameters [c6,c12] for each type pair, length n_types*n_types*2
 * * `n_types` - Number of atom types
 * * `box_data` - Box vectors (3 values for rectangular box, 9 for triclinic)
 * * `boundary_type` - 0=vacuum, 1=rectangular, 2=triclinic
 * * `crf_cut` - CRF cutoff distance
 * * `crf_2cut3i` - CRF parameter: 2/cutoff^3
 * * `crf_cut3i` - CRF parameter: 1/cutoff^3
 * * `forces` - Output: forces, flat array length n_atoms*3
 * * `energies` - Output: [e_lj, e_crf]
 * * `virial` - Output: virial tensor, flat array 3x3 = 9 elements
 */
void rust_lj_crf_innerloop(const float *positions,
                           const float *charges,
                           const unsigned int *iac,
                           unsigned int n_atoms,
                           const unsigned int *pairlist,
                           unsigned int n_pairs,
                           const double *lj_params,
                           unsigned int n_types,
                           const double *box_data,
                           unsigned int boundary_type,
                           double crf_cut,
                           double crf_2cut3i,
                           double crf_cut3i,
                           float *forces,
                           double *energies,
                           double *virial);

/**
 * C FFI: Parallel version of innerloop
 */
void rust_lj_crf_innerloop_parallel(const float *positions,
                                    const float *charges,
                                    const unsigned int *iac,
                                    unsigned int n_atoms,
                                    const unsigned int *pairlist,
                                    unsigned int n_pairs,
                                    const double *lj_params,
                                    unsigned int n_types,
                                    const double *box_data,
                                    unsigned int boundary_type,
                                    double crf_cut,
                                    double crf_2cut3i,
                                    double crf_cut3i,
                                    float *forces,
                                    double *energies,
                                    double *virial);

#endif /* GROMOS_RS_H */
