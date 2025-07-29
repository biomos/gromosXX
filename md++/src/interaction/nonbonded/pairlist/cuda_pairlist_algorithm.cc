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
 * @file cuda_pairlist_algorithm.cc
 * CUDA accelerated pairlist algorithm
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "math/periodicity.h"

#include "pairlist.h"
#include "pairlist_algorithm.h"
#include "cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

template<typename Backend>
interaction::CUDA_Pairlist_Algorithm<Backend>::CUDA_Pairlist_Algorithm()
: Pairlist_Algorithm(),
  m_solvent_solvent_timing(0.0) {}

/**
 * calculate center of geometries
 */
template<typename Backend>
int interaction::CUDA_Pairlist_Algorithm<Backend>::prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(0, "cuda pairlist algorithm : prepare");
  // m_impl.prepare(topo, conf, sim);
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  if (!sim.param().pairlist.atomic_cutoff) {

    // first put the chargegroups into the box
    m_impl.prepare_cog(conf, topo);

    // calculate cg cog's
    DEBUG(10, "calculating cg cog (" << topo.num_solute_chargegroups() << ")");
    m_cg_cog.resize(topo.num_solute_chargegroups());
    math::VArray const &pos = conf.current().pos;
    DEBUG(10, "pos.size() = " << pos.size());

    // calculate solute center of geometries
    topology::Chargegroup_Iterator
      cg1 =   topo.chargegroup_begin();
    
    unsigned int i = 0, num_cg = topo.num_solute_chargegroups();
    
    for(i=0; i < num_cg; ++cg1, ++i){
      cg1.cog(pos, m_cg_cog(i));
    }

  } // chargegroup based cutoff

  return 0;

}

template<typename Backend>
void interaction::CUDA_Pairlist_Algorithm<Backend>::update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       interaction::PairlistContainer & pairlist,
       unsigned int begin, unsigned int end,
       unsigned int stride)
{}

template<typename Backend>
void interaction::CUDA_Pairlist_Algorithm<Backend>::update_perturbed(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    interaction::PairlistContainer & pairlist,
    interaction::PairlistContainer & perturbed_pairlist,
    unsigned int begin, unsigned int end,
    unsigned int stride)
{}

// force instantiation
template class interaction::CUDA_Pairlist_Algorithm<util::cpuBackend>;

#ifdef USE_CUDA
template class interaction::CUDA_Pairlist_Algorithm<util::gpuBackend>;
#endif
