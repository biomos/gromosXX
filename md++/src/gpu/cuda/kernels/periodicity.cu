/**
 * @file hello_world.cu
 * @author poliak
 * example kernel function
 */

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"

#include "periodicity.h"

 __global__ void gpu::put_chargegroups_into_box(gpu::Topology gpu_topo,
                                                gpu::Configuration::View gpu_conf) {
//  __global__ void gpu::hello_world(float* a, float* b) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  // iterate over solute charge groups
  for (unsigned i = idx; i < gpu_topo.num_solute_chargegroups; i += blockDim.x * gridDim.x) {
    // iterate over atoms of the char
    const float3& my_pos = gpu_conf.current().pos[i];
  }
  // iterate over solvent charge groups
};

