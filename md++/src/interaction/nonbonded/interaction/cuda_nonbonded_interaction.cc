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
 * @file cuda_nonbonded_interaction.cc
 * template methods of CUDA_Nonbonded_Interaction.
 */
#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "simulation/parameter.h"

#include "interaction/interaction.h"
#include "interaction/interaction_types.h"
#include "interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "interaction/nonbonded/pairlist/pairlist.h"
#include "interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "interaction/nonbonded/interaction/storage.h"

#include "interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "interaction/nonbonded/interaction/nonbonded_set_interface.h"
#include "interaction/nonbonded/interaction/nonbonded_set.h"

#include "interaction/nonbonded/interaction/nonbonded_term.h"
#include "interaction/nonbonded/interaction/perturbed_nonbonded_term.h"

#include "interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"
#include "interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"

#include "interaction/nonbonded/interaction/perturbed_nonbonded_set.h"

#include "interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "interaction/nonbonded/interaction/cuda_nonbonded_interaction.h"

#include "util/debug.h"

#include "math/periodicity.h"
#include "math/boundary_checks.h"
#include "util/template_split.h"

#include "gpu/cuda/utils.h"

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::CUDA_Nonbonded_Interaction::CUDA_Nonbonded_Interaction(CUDA_Pairlist_Algorithm<util::gpuBackend> *pa)
: Nonbonded_Interaction(pa) {
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::CUDA_Nonbonded_Interaction::~CUDA_Nonbonded_Interaction() {
  delete m_pairlist_algorithm;
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::CUDA_Nonbonded_Interaction::
calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(4, "CUDA_Nonbonded_Interaction::calculate_interactions");

  m_timer.start(sim);

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled

  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;

  // std::cerr << "Nonbonded: steps = " << steps << std::endl;
  configuration::Configuration *p_conf = &conf;
  topology::Topology *p_topo = &topo;
  if ((sim.steps() % steps) == 0) {

    // std::cerr << "\tMULTISTEP: full non-bonded calculation" << std::endl;

    // shared memory do this only once
    if (m_pairlist_algorithm->prepare(*p_topo, *p_conf, sim))
      return 1;

    // have to do all from here (probably it's only one,
    // but then maybe it's clearer like it is...)
    // for (int i = 0; i < m_set_size; ++i) {
    //   if(m_nonbonded_set[i]->calculate_interactions(*p_topo, *p_conf, sim))
	  //     return 1;
    // }

    ///////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  } else {
    // std::cerr << "\tMULTISTEP: no non-bonded calculation" << std::endl;
  }

  DEBUG(6, "sets are done, adding things up...");
  store_set_data(*p_topo, *p_conf, sim);
  
  if (sim.param().multicell.multicell) {
    reduce_configuration(topo, conf, sim, *p_conf);
  }

  ////////////////////////////////////////////////////
  // printing pairlist
  ////////////////////////////////////////////////////
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))) {
    DEBUG(7, "print pairlist...");
    std::cerr << "printing pairlist!" << std::endl;
    print_pairlist(*p_topo, *p_conf, sim);
  }

  DEBUG(6, "CUDA_Nonbonded_Interaction::calculate_interactions done");
  m_timer.stop();

  return 0;
}

/**
 * initialize the arrays
 */
int interaction::CUDA_Nonbonded_Interaction::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  if (!quiet)
    os << "NONBONDED INTERACTION (CUDA)\n";

  // initialize data on GPU - positions, forces, energies
  // std::cuvector<float3> pos, force; // should be part of conf?
  // load all simulation and nonbonded parameters

  // pairlist then initializes pairlist
  

  configuration::Configuration * p_conf = &conf;
  topology::Topology * p_topo = &topo;

  if (!math::boundary_check_cutoff(p_conf->current().box, p_conf->boundary_type,
            sim.param().pairlist.cutoff_long)) {
    io::messages.add("box is too small: not twice the cutoff!",
            "configuration", io::message::error);
    return 1;
  }
  // if (m_pairlist_algorithm) {
  //     std::cout << "m_pairlist_algorithm address: " << m_pairlist_algorithm << std::endl;
  //     try {
  //         std::cout << "Runtime type: " << typeid(*m_pairlist_algorithm).name() << std::endl;
  //     } catch (const std::exception& e) {
  //         std::cerr << "Exception while accessing typeid: " << e.what() << std::endl;
  //     }
  // } else {
  //     std::cerr << "m_pairlist_algorithm is null!" << std::endl;
  // }

  // initialise the pairlist...
  m_pairlist_algorithm->init(*p_topo, *p_conf, sim, os, quiet);

  if (sim.param().nonbonded.method != simulation::el_reaction_field) {
    if (!quiet)
      os << "\tlattice-sum electrostatics\n";
    p_conf->lattice_sum().init(*p_topo, sim);
  }

  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.clear();

  // in case we do perturbation and eds at the same time, this is handled
  // in the eds_outer_loop. So we start with checking for EDS.
  if (sim.param().eds.eds || sim.param().perturbation.perturbation) {
    DEBUG(16, "creating EDS-perturbed nonbonded set");
    for (int i = 0; i < m_set_size; ++i) {
      m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
              m_parameter, i, m_set_size));
      DEBUG(16, "pushed back EDS-perturbed nonbonded set");
    }
  } else {
    for (int i = 0; i < m_set_size; ++i) {
      m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm,
              m_parameter, i, m_set_size));
    }
  }

  std::vector<Nonbonded_Set_Interface *>::iterator
  it = m_nonbonded_set.begin(),
          to = m_nonbonded_set.end();

  if (!quiet)
    os << "\tcreated " << m_nonbonded_set.size() << " set"
        <<  (m_nonbonded_set.size() > 1 ? "s" : "")    << "\n";

  bool q = quiet;
  for (; it != to; ++it) {
      (*it)->init(*p_topo, *p_conf, sim, os, q);
    // only print first time...
    q = true;
  }

  if (check_special_loop(*p_topo, *p_conf, sim, os, quiet) != 0) {
    io::messages.add("special solvent loop check failed", "CUDA_Nonbonded_Interaction",
            io::message::error);
    return 1;
  }
  if (!quiet)
    os << "END\n";
  DEBUG(9, "nonbonded init done");
  return 0;
}