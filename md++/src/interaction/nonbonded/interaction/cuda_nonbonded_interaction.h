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
 * @file cuda_nonbonded_interaction.h
 * GPU-accelerated LJ + reaction-field nonbonded interaction.
 *
 * Workflow per MD step:
 *   1. CPU Standard_Pairlist_Algorithm builds the pair list (handles
 *      exclusions, chargegroups, periodic images correctly).
 *   2. The CPU pair list is flattened to a CUDA-managed uint2 array.
 *   3. GPU zeros the force array and computes LJ-CRF forces in parallel.
 *   4. GPU forces are synchronised back to conf.current().force (CPU).
 *
 * The m_pairlist_algorithm member (CUDA_Pairlist_Algorithm<gpuBackend>)
 * is retained for future GPU pairlist acceleration; it is not used in
 * this first vanilla implementation.
 */

#pragma once

#include "interaction.h"
#include "nonbonded_parameter.h"
#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "gpu/cuda/memory/cuvector.h"

namespace topology   { class Topology; }
namespace configuration { class Configuration; }
namespace simulation  { class Simulation; }

namespace interaction
{

/**
 * @class CUDA_Nonbonded_Interaction
 * Computes LJ + reaction-field nonbonded forces on the GPU.
 */
class CUDA_Nonbonded_Interaction : public Nonbonded_Interaction
{
public:

    /**
     * Constructor.
     * @param pa  GPU pairlist algorithm (reserved for future GPU pairlist).
     */
    explicit CUDA_Nonbonded_Interaction(
        CUDA_Pairlist_Algorithm<util::gpuBackend>* pa);

    virtual ~CUDA_Nonbonded_Interaction();

    /** Calculate nonbonded interactions using GPU force kernels. */
    virtual int calculate_interactions(topology::Topology&        topo,
                                       configuration::Configuration& conf,
                                       simulation::Simulation&    sim);

    /** Initialise GPU data structures (LJ matrix, CRF constants, pairlist). */
    virtual int init(topology::Topology&        topo,
                     configuration::Configuration& conf,
                     simulation::Simulation&    sim,
                     std::ostream&              os    = std::cout,
                     bool                       quiet = false);

private:
    /** CPU pairlist algorithm – builds pairs with correct exclusions. */
    Standard_Pairlist_Algorithm m_std_pairlist_alg;

    /** CPU-side pairlist container populated each step (or every skip_step). */
    PairlistContainer m_pairlist;

    /** GPU LJ parameter matrix (indexed by integer atom code). */
    gpu::LJParams m_lj_params;

    /** CRF constants and cutoffs for GPU force kernels. */
    gpu::NbSimParams m_nb_params;

    /** Flat GPU pair list; each element is {atom_i, atom_j}. */
    gpu::cuvector<uint2> m_gpu_pairs;

    /** Convert CPU pairlist into m_gpu_pairs. */
    void build_flat_pairlist(const PairlistContainer& pl);
};

} // namespace interaction
