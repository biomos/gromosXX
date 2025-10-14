


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

#include "gpu/cuda/utils.h"

#define NUM_THREADS_PER_BLOCK 256

interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::CUDA_Pairlist_Algorithm_Impl()
        : CUDA_Pairlist_Algorithm_ImplBase() {
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

void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::prepare_cog(
                                configuration::Configuration & conf,
                                topology::Topology & topo,
                                simulation::Simulation & sim) {
    set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
    
    if (!sim.param().pairlist.atomic_cutoff){
        conf.copy_to_gpu();
        const size_t num_solute_cg = topo.num_solute_chargegroups();
        const size_t num_cg = topo.num_chargegroups();
        m_cg_cog.resize(num_solute_cg);
        m_cg_cells.resize(num_cg);
        SPLIT_BOUNDARY(_prepare_cog, conf, topo);
    }
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::_prepare_cog<B>(
                                    configuration::Configuration & conf,
                                    topology::Topology & topo) {
    DEBUG(10, "putting chargegroups into box");

    dim3 dimGrid(1);
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();
    gpu::Periodicity<B> periodicity(conf.current().box);
    std::cout << "m_cutoff_long: " << m_cutoff_long << std::endl;
    periodicity.set_cell_size(m_cutoff_long);
    gpu::prepare_cog_kernel<<<dimGrid, dimBlock>>>(topo.get_gpu_view(),
                                                    conf.get_gpu_view(),
                                                    periodicity,
                                                    m_cg_cog.view(),
                                                    m_cg_cells.view());
};

void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::set_cutoff(
                    double const cutoff_short,
                    double const cutoff_long)
{
    m_cutoff_long = cutoff_long;
    m_cutoff_short = cutoff_short;
    m_cutoff_short_2 = cutoff_short * cutoff_short;
    m_cutoff_long_2  = cutoff_long * cutoff_long;
};

void interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>::reorder(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{
    
};

// explicit instantiation
template class interaction::CUDA_Pairlist_Algorithm_Impl<util::gpuBackend>;
