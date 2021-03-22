/**
 * @file gathering.cc gathering for QM/MM simulations
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../math/periodicity.h"
#include "../../../util/template_split.h"

// special interactions
#include "../../../interaction/special/qmmm/mm_atom.h"

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
  // new: ref = qm_pos(0)
  const math::Vec & ref = qm_pos(0);
  math::Vec r;
  for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
	// gather the MM atoms with respect to qm_pos(0)
 	periodicity.nearest_image(ref, mm_atoms[i].pos, r);
     mm_atoms[i].pos = ref - r;

    
  }
}

void interaction::gather_mm_atoms(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        math::VArray & qm_pos,
        std::vector<interaction::MM_Atom>& mm_atoms) {
  SPLIT_BOUNDARY(_gather_mm_atoms, topo, conf, sim, qm_pos, mm_atoms);
}
