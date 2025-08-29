/**
 * @file hello_world.cu
 * @author poliak
 * example kernel function
 */

#include "gpu/cuda/cuheader.h"


#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"

#include "gpu/cuda/math/periodicity.h"

#include "periodicity.h"

 __global__ void gpu::put_chargegroups_into_box(gpu::Topology gpu_topo,
                                                gpu::Configuration::View gpu_conf) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  // iterate over solute charge groups
  for (unsigned i = idx; i < gpu_topo.num_solute_chargegroups; i += blockDim.x * gridDim.x) {
    // iterate over atoms of the char
    const float3& my_pos = gpu_conf.current().pos[i];
  }
  // iterate over solvent charge groups
};

template <typename VecType>
__device__ __forceinline__ VecType gpu::calculate_cog(VecType* pos, int begin, int end) {
  VecType cog{0.,0.,0.};
    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  for (int a_i = begin; a_i < end; ++a_i) {
    printf("idx: %d, a_i: %d = (%f,%f,%f)\n", idx, a_i, pos[a_i].x, pos[a_i].y, pos[a_i].z);
    cog += pos[a_i];
  }
  return cog;
}

template <math::boundary_enum BOUNDARY>
__global__ void gpu::prepare_cog_kernel(gpu::Topology topo,
                                        gpu::Configuration::View conf,
                                        gpu::Periodicity<BOUNDARY> periodicity) {
    const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned stride = blockDim.x * gridDim.x;
    unsigned cg_i = idx;
    for(; cg_i < topo.num_chargegroups; cg_i += stride) {
        const int cg_begin = topo.chargegroup[cg_i];
        const int cg_end = topo.chargegroup[cg_i+1];
        FPL3_TYPE* pos = conf.current().pos;

        FPL3_TYPE cog;
        if (cg_i < topo.num_solute_chargegroups) {
          cog = calculate_cog(pos, cg_begin, cg_end);
        } else { // solvent
          cog = pos[cg_begin];
        }
        printf("cog cg: %d = (%f,%f,%f)\n", cg_i, cog.x, cog.y, cog.z);
        FPL3_TYPE v_box = cog;
        periodicity.put_into_box(v_box);
        FPL3_TYPE trans = v_box - cog;
        for (int a_i = cg_begin; a_i < cg_end; ++a_i) {
            assert(conf.current().size > a_i);
            conf.current().pos[a_i] += trans;
        } // loop over atoms
    } // loop over all cg's
};

template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology topo, gpu::Configuration::View conf, gpu::Periodicity<math::vacuum> periodicity);

template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology topo, gpu::Configuration::View conf, gpu::Periodicity<math::rectangular> periodicity);

template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology topo, gpu::Configuration::View conf, gpu::Periodicity<math::triclinic> periodicity);


