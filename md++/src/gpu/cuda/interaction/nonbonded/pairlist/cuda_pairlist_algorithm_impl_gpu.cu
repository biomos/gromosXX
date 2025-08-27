


#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

// #include "math/periodicity.h"
#include "math/boundary_implementation.h"
#include "gpu/cuda/math/periodicity_gpu.h"

#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/kernels/hello_world.h"

#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl.h"
#include "cuda_pairlist_algorithm_impl_gpu.h"

#include "util/debug.h"

#define NUM_THREADS_PER_BLOCK 256

interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::CUDA_Pairlist_Algorithm_Impl() {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl<util::gpuBackend> constructor");
};

int interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::init(topology::Topology &topo, 
    configuration::Configuration &conf,
    simulation::Simulation &sim,
    std::ostream &os,
    bool quiet) {
    // GPU-specific code
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::init");
    return 0;
};

void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::
                    prepare_cog(configuration::Configuration & conf,
                                topology::Topology & topo) {
    SPLIT_BOUNDARY(_prepare_cog, conf, topo);
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::
                    _prepare_cog<B>(configuration::Configuration & conf,
                                    topology::Topology & topo) {
    DEBUG(10, "putting chargegroups into box");
    gpu::PeriodicityGpu<B> periodicity(conf.current().box);
    // periodicity.put_chargegroups_into_box(conf, topo);
    periodicity.put_chargegroups_into_box(conf, topo);

    dim3 dimGrid(1);
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();

    // DEBUG(0, "topo.get_gpu_view().num_atoms: " << topo.get_gpu_view().num_atoms);
    // cudaMemcpy(d_a, arr, sizeof(float)*2, cudaMemcpyHostToDevice);
    gpu::hello_world<<<dimGrid, dimBlock>>>(topo.get_gpu_view(), conf.get_gpu_view());
    // gpu::prepare_cog<<<dimGrid, dimBlock>>>(topo.get_gpu_view(), conf.get_gpu_view());
    // gpu::hello_world<<<dimGrid, dimBlock>>>(d_a,d_b);
    cudaDeviceSynchronize();
    // cudaMemcpy(arr, d_a, sizeof(float)*2, cudaMemcpyDeviceToHost);
    // DEBUG(0, "a: " << a << ", b: " << b);
    // cudaFree(d_a);
}

// explicit instantiation
template class interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>;
