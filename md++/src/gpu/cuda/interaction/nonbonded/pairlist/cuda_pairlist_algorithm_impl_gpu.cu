


#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

// #include "math/periodicity.h"
#include "math/boundary_implementation.h"
#include "gpu/cuda/math/periodicity.h"

#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/kernels/periodicity.h"

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

    conf.copy_to_gpu();
    SPLIT_BOUNDARY(_prepare_cog, conf, topo);
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::
                    _prepare_cog<B>(configuration::Configuration & conf,
                                    topology::Topology & topo) {
    DEBUG(10, "putting chargegroups into box");


    dim3 dimGrid(1);
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();
    gpu::Periodicity<B> periodicity(conf.current().box);
    gpu::prepare_cog_kernel<<<dimGrid, dimBlock>>>(topo.get_gpu_view(), conf.get_gpu_view(), periodicity);

    cudaDeviceSynchronize();
};

// explicit instantiation
template class interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>;
