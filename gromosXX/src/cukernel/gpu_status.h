/**
 * @file gpu_status.h
 * Defines all the pointers and other importent info for communication between
 * host (CPU) and device (GPU)
 */

#ifndef _GPU_STATUS_H
#define	_GPU_STATUS_H

#include "parameter.h"
#include "pairlist.h"
#include "lib/types.h"

/**
 * @struct gpu_status
 * Defines all the pointers and other importent info for communication between
 * host (CPU) and device (GPU)
 */
struct gpu_status {
  int device;
  cudakernel::simulation_parameter host_parameter;
  float2 * host_energy;
  float3 * host_forces;
  float9 * host_virial;
  float3 * host_pos;
  float3 * host_new_pos;
  float3 * host_old_pos;
  double3 * host_double_new_pos;
  double3 * host_double_old_pos;

  float3 * dev_pos;
  float3 * dev_old_pos;
  float3 * dev_new_pos;
  double3 * dev_double_new_pos;
  double3 * dev_double_old_pos;
  float3 * dev_forces;
  float9 * dev_virial;
  float9 * dev_factor;
  double9 * dev_double_factor;
  cudakernel::pairlist dev_pl_short, dev_pl_long;
  cudakernel::simulation_parameter * dev_parameter;
  cudakernel::lj_crf_parameter * dev_lj_crf_parameter;
  float2 * dev_energy;
  float3 * dev_const_length2;
  double3 * dev_double_const_length2;
  float * dev_tol;
  double * dev_double_tol;
  float3 * dev_mass;
  double3 * dev_double_mass;
  int * dev_shake_fail_mol;
  unsigned int * dev_highest_index;

  cudakernel::constraint * dev_constr;
};

#endif	/* _GPU_STATUS_H */

