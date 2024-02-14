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

