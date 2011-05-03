/**
 * @file local_elevation_interaction.cc
 * apply LEUS
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/special/qmmm/qm_storage.h>
#include <interaction/special/qmmm/qm_worker.h>
#include <interaction/special/qmmm/mndo_worker.h>
#include <interaction/special/qmmm_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>
#include <vector>
#include <map>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

interaction::QMMM_Interaction::~QMMM_Interaction() {
  if (worker != NULL)
    delete worker;
}

int interaction::QMMM_Interaction::
calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  m_timer.start();
  storage.zero();
  
  // this might fail due to user input!
  if (worker->run_QM(topo, conf, sim, storage))
    return 1;

  // add QM forces
  for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
    conf.current().force(i) += storage.force(i);
  }
  
  conf.current().energies.qm_total = storage.energy;
  m_timer.stop();
 
  return 0;
}

int interaction::QMMM_Interaction::init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os,
            bool quiet) {
  storage.resize(topo.num_atoms(), topo.qm_zone().size());
  
  worker = interaction::QM_Worker::get_instance(sim);
  if (worker == NULL) return 1;
  if (worker->init(topo, conf, sim)) return 1;

  if (!quiet) {
    os << "QM/MM" << std::endl
            << "  QM worker: " << worker->name() << std::endl;

    os << "  QM zone: " << std::endl;
    for (std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
      os << "     - " << std::setw(5) << it->index+1 << std::setw(5) << it->atomic_number;
      if (it->link) {
        os << " link atom";
      }
      os << std::endl;
    }
    os << "END" << std::endl;
  }
  
  // set charge of QM atoms to 0
  for(std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
    topo.charge()[it->index] = 0.0;
  }
  
  // remove bonded terms
  for (unsigned int i = 0; i < topo.solute().bonds().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().bonds()[i].i) || 
        topo.in_qm_zone(topo.solute().bonds()[i].j)) {
      topo.solute().bonds().erase(topo.solute().bonds().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().cgbonds().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().cgbonds()[i].i) || 
        topo.in_qm_zone(topo.solute().cgbonds()[i].j)) {
      topo.solute().cgbonds().erase(topo.solute().cgbonds().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().distance_constraints().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().distance_constraints()[i].i) || 
        topo.in_qm_zone(topo.solute().distance_constraints()[i].j)) {
      topo.solute().distance_constraints().erase(topo.solute().distance_constraints().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().angles().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().angles()[i].i) || 
        topo.in_qm_zone(topo.solute().angles()[i].j) || 
        topo.in_qm_zone(topo.solute().angles()[i].k)) {
      topo.solute().angles().erase(topo.solute().angles().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().dihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().dihedrals()[i].i) || 
        topo.in_qm_zone(topo.solute().dihedrals()[i].j) || 
        topo.in_qm_zone(topo.solute().dihedrals()[i].k) ||
        topo.in_qm_zone(topo.solute().dihedrals()[i].l)) {
      topo.solute().dihedrals().erase(topo.solute().dihedrals().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().improper_dihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().improper_dihedrals()[i].i) || 
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].j) || 
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].k) ||
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].l)) {
      topo.solute().improper_dihedrals().erase(topo.solute().improper_dihedrals().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().crossdihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().crossdihedrals()[i].a) || 
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].b) || 
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].c) ||
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].d) ||
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].e) ||
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].f) ||
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].g) ||
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].h)) {
      topo.solute().crossdihedrals().erase(topo.solute().crossdihedrals().begin() + i);
      --i;
    }
  }
  return 0;
}

