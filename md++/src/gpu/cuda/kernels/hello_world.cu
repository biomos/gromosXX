/**
 * @file hello_world.cu
 * @author poliak
 * example kernel function
 */

// #include "../cuheader.h"
// #include "../memory/topology_struct.h"
// #include "../memory/configuration_state.h"

#include "hello_world.h"

//  __global__ void gpu::hello_world(gpu::topology_struct& topo_struct, gpu::configuration_state_struct& conf_state_struct) {
 __global__ void gpu::hello_world(float* a, float* b) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  // for (unsigned i = idx; i < topo_struct.num_atoms; i += blockDim.x * gridDim.x) {
  //   Perform some operation on the configuration state
  //   float3* my_pos = &conf_state_struct.pos[i]
  //   printf("threadIdx %d, blockIdx %d: processing %d: pos: (%f,%f,%f) iac: %d\n", threadIdx.x, blockIdx.x, i, my_pos->x, my_pos->y, my_pos->z, topo.iac[i]);
  atomicAdd(a, *b);
  atomicAdd(b, *a);
  // }
};

