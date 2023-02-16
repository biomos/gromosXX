/**
 * @file cudaShake.h
 * m_shake algorithm
 */

#ifndef CUKERNEL_CONSTRAINTS
#include "parameter.h"
#include "macros.h"

namespace cudakernel {
  /**
   * solve the constraints using the M-SHAKE algorithm
   * @param[inout] new_pos the new positions
   * @param[in] old_pos the old positions
   * @param[in] dev_params the simulation parameters
   * @param[out] shake_file_mol. The molecule for which the constraints algorithm failed. -1 for success.
   * @param[in] tol the tolerance
   * @param[in] mass the masses fo the atoms in a molecule
   * @param[in] const_length2 the constraint lengths squared
   * @param[in] factor the constraint matrix
   * @param[in] highest_mol_index the hightest index of the molecule in the array
   */
  __global__ void kernel_Calc_Shake
  (
          VECTOR * new_pos, VECTOR * old_pos,
          cudakernel::simulation_parameter * dev_params,
          int *shake_fail_mol, FL_PT_NUM * tol, 
          VECTOR * mass,
          VECTOR * const_length2,
          MATRIX * factor, unsigned int highest_mol_index
          );

}

#endif

