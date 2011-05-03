/**
 * @file gathering.cc gathering for QM/MM simulations
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>
#include <util/template_split.h>

// special interactions
#include <interaction/special/qmmm/mm_atom.h>

#include "gathering.h"

template<math::boundary_enum b>
void _gather_qmzone(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        math::VArray & qm_pos) {
  math::Periodicity<b> periodicity(conf.current().box);
  const math::Vec & ref = qm_pos(0);
  math::Vec r;
  for (unsigned int i = 1; i < topo.qm_zone().size(); ++i) {
    periodicity.nearest_image(ref, qm_pos(i), r);
    qm_pos(i) = ref - r;
  }
}

void interaction::gather_qmzone(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        math::VArray & qm_pos) {
  SPLIT_BOUNDARY(_gather_qmzone, topo, conf, sim, qm_pos);
}

template<math::boundary_enum b>
void _gather_mm_atoms(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        math::VArray & qm_pos,
        std::vector<interaction::MM_Atom>& mm_atoms) {
  math::Periodicity<b> periodicity(conf.current().box);
  
  math::Vec r;
  for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
    // determine the minimal distance to all QM atoms
    // and use the closest as the reference for gathering
    math::Vec min_ref(qm_pos(0));
    periodicity.nearest_image(qm_pos(0), mm_atoms[i].pos, r);
    math::Vec min_r(r);
    double min_dist = math::abs2(r);
    for(unsigned int j = 1; j < qm_pos.size(); ++j) {
      periodicity.nearest_image(qm_pos(j), mm_atoms[i].pos, r);
      double r2 = math::abs2(r);
      if (r2 < min_dist) {
        min_ref = qm_pos(j);
        min_dist = r2;
        min_r = r;
      }
    }
    
    // gather the atom
    mm_atoms[i].pos = min_ref - min_r;
  }
}

void interaction::gather_mm_atoms(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        math::VArray & qm_pos,
        std::vector<interaction::MM_Atom>& mm_atoms) {
  SPLIT_BOUNDARY(_gather_mm_atoms, topo, conf, sim, qm_pos, mm_atoms);
}