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

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "../../../util/debug.h"
#include "../../../util/template_split.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::CUDA_Pairlist_Algorithm::CUDA_Pairlist_Algorithm()
: interaction::Pairlist_Algorithm(),
  m_solvent_solvent_timing(0.0) {}

/**
 * calculate center of geometries
 */
int interaction::CUDA_Pairlist_Algorithm::prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  return 0;
}

void interaction::CUDA_Pairlist_Algorithm::update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       interaction::PairlistContainer & pairlist,
       unsigned int begin, unsigned int end,
       unsigned int stride)
{}

void interaction::CUDA_Pairlist_Algorithm::update_perturbed(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    interaction::PairlistContainer & pairlist,
    interaction::PairlistContainer & perturbed_pairlist,
    unsigned int begin, unsigned int end,
    unsigned int stride)
{}