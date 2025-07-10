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
 * @file parameter.h
 * the device parameters copied preferably to the constant memory
 */
#ifndef INCLUDED_CUKERNEL_PARAMETER_H
#define INCLUDED_CUKERNEL_PARAMETER_H

#define MAX_ATOMS_SOLVENT 6

namespace cukernel {
  /**
   * @struct lj_crf_parameter
   * struct containing the interaction parameters for pairwise forces
   * this struct is used to build a NxN-interaction parameter matrix
   * where N is the number of different atom-types
   */
  struct lj_crf_parameter {
    /**
     * the Lennard-Jones C6 constant of the pair
     */
    float c6;
    /**
     * the Lennard-Jones C12 constant of the pair
     */
    float c12;
    /**
     * product of charges times 1/(4 pi eps0)
     */
    float q;
  };

  /**
   * @struct simulation_parameter
   *
   * struct containing simulation parameters needed by kernels
   */
  struct simulation_parameter {
    static const unsigned max_atoms_solvent = MAX_ATOMS_SOLVENT;
    /**
     * precalculated lj and crf products
     */
    lj_crf_parameter solvent_lj_crf[max_atoms_solvent*max_atoms_solvent];

    struct box_struct {
      /**
       * the box edges
       */
      float3 full;
      /**
       * the inverted box edges
       */
      float3 inv;
      /**
       * half the box edges
       */
      float3 half;
    } box;
    struct num_atoms_struct {
      /**
       * number of all atoms
       */
      unsigned total;
      /**
       * number of solute atoms
       */
      unsigned solute;
      /**
       * number of solvent atoms
       */
      unsigned solvent;
    } num_atoms;
    /**
     * the long-range cutoff
     */
    float cutoff_long;
    /**
     * the squared long-range cutoff
     */
    float cutoff_long_2;
    /**
     * the short-range cutoff
     */
    float cutoff_short;
    /**
     * the squared short-range cutoff
     */
    float cutoff_short_2;
    /**
     * reaction field constant
     */
    float crf_2cut3i;
    /**
     * reaction field constant
     */
    float crf_cut;
    /**
     * reaction field constant
     */
    float crf_cut3i;
    /**
     * the number of atoms per solvent molecule
     */
    unsigned int num_atoms_per_mol;
    /**
     * the number of solvent molecules
     */
    unsigned int num_solvent_mol;
    /**
     * the estimated size of the long-range neighbor list
     */
    unsigned int estimated_neighbors_long;
    /**
     * the estimated size of the short-range neighbor list
     */
    unsigned int estimated_neighbors_short;
    /**
     * The number of GPUs
     */
    unsigned int num_of_gpus;
    /**
     * The current number of the gpus
     */
    unsigned int gpu_id;
  };

  /**
   * @struct constraint
   * Used to copy all the informations from the two_body_term_struct and from
   * paramter to the GPU
   */
  struct constraint {

    constraint() : i(0), j(0), length(0.0) {
    }

    constraint(unsigned int i, unsigned int j, double length) : i(i), j(j), length(length) {
    }

    void set_values(unsigned int ip, unsigned int jp, double r0) {
      i = ip;
      j = jp;
      length = r0;
    }
    unsigned int i, j;
    double length;
  };
}
#endif

