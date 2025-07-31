/**
 * @file hello_world.cu
 * @author poliak
 * example kernel function
 */

#include "gpu/cuda/cuheader.h"
#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"

#include "hello_world.h"

 __global__ void gpu::hello_world(gpu::topology_struct topo, gpu::Configuration::View conf) {
//  __global__ void gpu::hello_world(float* a, float* b) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  for (unsigned i = idx; i < topo.num_atoms; i += blockDim.x * gridDim.x) {
    // Perform some operation on the configuration state
    const float3& my_pos = conf.current().pos[i];
    printf("threadIdx %d, blockIdx %d: processing %d: pos: (%f,%f,%f) iac: %d\n", threadIdx.x, blockIdx.x, i, my_pos.x, my_pos.y, my_pos.z, topo.iac[i]);
  }
};

