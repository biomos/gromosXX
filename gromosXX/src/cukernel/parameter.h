/**
 * @file parameter.h
 * the accelerator's parameters
 */
#ifndef INCLUDED_CUKERNEL_PARAMETER_H
#define INCLUDED_CUKERNEL_PARAMETER_H

namespace cudakernel {

  /**
   * @struct simulation_parameter
   *
   * struct containing simulation parameters needed by kernels
   * initialized in cudaInit(...)
   */
  struct simulation_parameter {
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
     * the box edges
     */
    float box_x;
    float box_y;
    float box_z;
    /**
     * the inverted box edges
     */
    float box_inv_x;
    float box_inv_y;
    float box_inv_z;
    /**
     * half the box edges
     */
    float box_half_x;
    float box_half_y;
    float box_half_z;
    
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
     * number of solvent atoms
     */
    unsigned int num_atoms;
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
     * The current number of the gpu
     */
    unsigned int gpu_id;
  };

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

