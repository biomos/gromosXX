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

template <math::boundary_enum BOUNDARY>
__global__ void gpu::put_chargegroups_into_box_kernel(
    gpu::Topology::View topo,
    gpu::Configuration::View conf,
    gpu::Periodicity<BOUNDARY> periodicity) {

    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned stride = blockDim.x * gridDim.x;

    for (unsigned cg_i = idx; cg_i < topo.num_chargegroups; cg_i += stride) {
        unsigned cg_begin = topo.chargegroup[cg_i];
        unsigned cg_end   = topo.chargegroup[cg_i+1];
        periodicity.prepare_chargegroup(
            cg_i, cg_begin, cg_end,
            conf.current().pos,
            topo.num_solute_chargegroups,
            nullptr,  // no cg_cog
            nullptr   // no cg_cells
        );
    }
}

template <math::boundary_enum BOUNDARY>
__global__ void gpu::prepare_cog_kernel(
    gpu::Topology::View topo,
    gpu::Configuration::View conf,
    gpu::Periodicity<BOUNDARY> periodicity,
    math::CuVArray::View cg_cog,
    gpu::cuvector<ushort4>::View cg_cells) {

    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned stride = blockDim.x * gridDim.x;

    for (unsigned cg_i = idx; cg_i < topo.num_chargegroups; cg_i += stride) {
        unsigned cg_begin = topo.chargegroup[cg_i];
        unsigned cg_end   = topo.chargegroup[cg_i+1];
        periodicity.prepare_chargegroup(
            cg_i, cg_begin, cg_end,
            conf.current().pos,
            topo.num_solute_chargegroups,
            &cg_cog,   // store cog
            &cg_cells // store cell
        );
    }
}

// explicit instantiation to allow linking
template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology::View topo,
    gpu::Configuration::View conf,
    gpu::Periodicity<math::vacuum> periodicity,
    math::CuVArray::View cg_cog,
    gpu::cuvector<ushort4>::View cg_cells);

template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology::View topo,
    gpu::Configuration::View conf,
    gpu::Periodicity<math::rectangular> periodicity,
    math::CuVArray::View cg_cog,
    gpu::cuvector<ushort4>::View cg_cells);

template __global__ void gpu::prepare_cog_kernel(
    gpu::Topology::View topo,
    gpu::Configuration::View conf,
    gpu::Periodicity<math::triclinic> periodicity,
    math::CuVArray::View cg_cog,
    gpu::cuvector<ushort4>::View cg_cells);


