/**
 * @file local_elevation_interaction.cc
 * apply LEUS
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/special/qmmm/qm_storage.h"
#include "../../interaction/special/qmmm/mm_atom.h"
#include "../../interaction/special/qmmm/qm_worker.h"
#include "../../interaction/special/qmmm/mndo_worker.h"
#include "../../interaction/special/qmmm_interaction.h"
#include "../../interaction/special/qmmm/gathering.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"
#include "../../math/boundary_checks.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

interaction::QMMM_Interaction::~QMMM_Interaction() {
  if (worker != NULL)
    delete worker;
}

int interaction::QMMM_Interaction::prepare(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  // determine which MM atoms to include as point charges
  m_timer.start("MM atoms determination");
  interaction::determine_mm_atoms(topo, conf, sim, mm_atoms);
  m_timer.stop("MM atoms determination");

  // get the position of the QM atoms and gather
  m_timer.start("gathering");
  {
    qm_pos.resize(topo.qm_zone().size());
    unsigned int i = 0;
    for (std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++i) {
      qm_pos(i) = conf.current().pos(it->index);
    }
  }
  interaction::gather_qmzone(topo, conf, sim, qm_pos);
  interaction::gather_mm_atoms(topo, conf, sim, qm_pos, mm_atoms);
  m_timer.stop("gathering");
  
  return 0;
}

int interaction::QMMM_Interaction::
calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  m_timer.start();
  storage.zero();
  
  // if polarisation is in use this was already called.
  if (!sim.param().polarise.cos) 
    prepare(topo, conf, sim);

  // this might fail due to user input!
  m_timer.start("QM worker");
  if (worker->run_QM(topo, conf, sim, qm_pos, mm_atoms, storage)) {
    m_timer.stop("QM worker");
    return 1;
  }
  m_timer.stop("QM worker");

  // add QM forces
  for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
    conf.current().force(i) += storage.force(i);
    if (topo.is_polarisable(i)) {
      conf.current().force(i) += storage.cos_force(i);
    }
  }
  
  conf.current().energies.qm_total = storage.energy;
  m_timer.stop();
 
  return 0;
}

int interaction::QMMM_Interaction::add_electric_field_contribution(topology::Topology& topo,
        configuration::Configuration& conf, simulation::Simulation& sim, 
        math::VArray & electric_field) {
  
  m_timer.start();
  m_timer.start("polarisation");
  
  storage.zero();
  if (worker->run_QM(topo, conf, sim, qm_pos, mm_atoms, storage)) {
    return 1;
  }
  for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
    if (!topo.is_polarisable(i))
      continue;

    // get electric field at either the charge or the COS site.
    math::Vec e;
    switch (sim.param().polarise.efield_site) {
      case simulation::ef_atom:
        e = storage.force(i) / (topo.charge(i) - topo.coscharge(i));
        break;
      case simulation::ef_cos:
        e = storage.cos_force(i) / topo.coscharge(i);
        break;
      default:
        io::messages.add("Electric field site not implemented.",
                "QMMM_Interaction", io::message::critical);
        return 1;
    }
    // add the electric field contribution
    electric_field(i) += e;
  }
  
  m_timer.stop("polarisation");
  m_timer.stop();
  return 0;
}

int interaction::QMMM_Interaction::init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os,
            bool quiet) { 
  if (topo.qm_zone().empty()) {
    io::messages.add("No QM zone defined", "QMMM_Interaction", io::message::error);
    return 1;
  }
  
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
    os.precision(3);
    if (sim.param().qmmm.cutoff > 0.0) {
      os << "  Cutoff for inclusion of MM atoms: " << std::setw(8)
              << sim.param().qmmm.cutoff << std::endl;
    } else {
      os << "  All atoms are included as MM atoms." << std::endl;
    }
    os << "END" << std::endl;
  }
  
  if (sim.param().qmmm.cutoff > 0.0 && 
          !math::boundary_check_cutoff(conf.current().box, conf.boundary_type,
          sim.param().qmmm.cutoff)) {
    io::messages.add("The cutoff RCUTQ is too large for the size of the "
            "computational box.", "QMMM_Interaction", io::message::error);
    return 1;
  }
  
  if (sim.param().qmmm.cutoff == 0.0 && conf.boundary_type != math::vacuum) {
    io::messages.add("A cutoff RCUTQ >= 0.0 has to be used for simulations "
            "using periodic boundary conditions.", "QMMM_Interaction", 
            io::message::error);
    return 1;
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

