/**
 * @file cuda_pairlist_algorithm_impl.cu
 * GPU-native implementation of CUDA_Pairlist_Algorithm_Impl.
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "math/boundary_implementation.h"
#include "gpu/cuda/math/periodicity.h"

#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/kernels/periodicity.h"

#include "cuda_pairlist_algorithm_impl.h"

#include "gpu/cuda/utils.h"

#define NUM_THREADS_PER_BLOCK 256

interaction::CUDA_Pairlist_Algorithm_Impl::CUDA_Pairlist_Algorithm_Impl() {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl constructor");
};

int interaction::CUDA_Pairlist_Algorithm_Impl::init(topology::Topology &topo,
    configuration::Configuration &conf,
    simulation::Simulation &sim,
    std::ostream &os,
    bool quiet) {
    DEBUG(0, "CUDA_Pairlist_Algorithm_Impl::init");
    return 0;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::set_cutoff(
                    double const cutoff_short,
                    double const cutoff_long)
{
    m_cutoff_long = cutoff_long;
    m_cutoff_short = cutoff_short;
    m_cutoff_short_2 = cutoff_short * cutoff_short;
    m_cutoff_long_2  = cutoff_long * cutoff_long;
};

void interaction::CUDA_Pairlist_Algorithm_Impl::prepare_cog(
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
        m_cg_sort_key.resize(num_cg);
        SPLIT_BOUNDARY(_prepare_cog, conf, topo);
    }
}

template<math::boundary_enum B>
void interaction::CUDA_Pairlist_Algorithm_Impl::_prepare_cog<B>(
                                    configuration::Configuration & conf,
                                    topology::Topology & topo) {
    DEBUG(10, "putting chargegroups into box");

    const unsigned num_cg = static_cast<unsigned>(topo.num_chargegroups());
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((num_cg + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();
    gpu::Periodicity<B> periodicity(conf.current().box);
    periodicity.set_cell_size(m_cutoff_long);
    gpu::prepare_cog_kernel<<<dimGrid, dimBlock>>>(topo.get_gpu_view(),
                                                    conf.get_gpu_view(),
                                                    periodicity,
                                                    m_cg_cog.view(),
                                                    m_cg_cells.view(),
                                                    m_cg_sort_key.view());
};

void interaction::CUDA_Pairlist_Algorithm_Impl::reorder(
                configuration::Configuration & conf,
                topology::Topology & topo,
                simulation::Simulation & sim)
{

};
