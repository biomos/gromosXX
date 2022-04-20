/**
 * @file xtb_worker.cc
 * interface to the XTB software package
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"

#include "../../../io/blockinput.h"

#include "../../../util/timing.h"
#include "../../../util/system_call.h"
#include "../../../util/debug.h"

#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "xtb_worker.h"

#include "xtb.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::XTB_Worker::XTB_Worker() : QM_Worker("XTB Worker"), param(nullptr) {}; 

interaction::XTB_Worker::~XTB_Worker() {
  // release Fortran resources
  xtb_delMolecule(&mol);
  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delEnvironment(&env);
}

int interaction::XTB_Worker::init(const topology::Topology& topo
                                , const configuration::Configuration& conf
                                , simulation::Simulation& sim
                                , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());

  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.xtb);
  QM_Worker::param = this->param;

  for (std::set<MM_Atom>::const_iterator 
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    std::cout << it->index << std::endl;
    std::cout << it->atomic_number << std::endl;
    std::cout << it->charge << std::endl;
    std::cout << math::v2s(it->pos) << std::endl;
  }

  std::cout << "size iac() : " << topo.iac().size() << std::endl;
  for (unsigned int i = 0; i < topo.iac().size(); ++i) {
    std::cout << "atom with index: " << i << " has IAC code: " << topo.iac(i) + 1 << " " << topo.is_qm(i) << std::endl;
  }

  std::cout << "Size of qm_atomic_number: " << topo.qm_atomic_number().size() << std::endl;
  int items_not_zero = std::count_if(topo.qm_atomic_number().begin(), topo.qm_atomic_number().end(), [](int item) {return item != 0;});
  int items_zero = topo.qm_atomic_number().size() - items_not_zero;
  std::cout << "items not zero: " << items_not_zero << std::endl;
  std::cout << "items zero: " << items_zero << std::endl;
  for (const auto& e : topo.qm_atomic_number()) {
    std::cout << e << std::endl;
  }

  this->charge = qm_zone.charge() * this->param->unit_factor_charge;
  this->uhf = (qm_zone.spin_mult() - 1) / 2; // spin multiplicity to # unpaired electrons

  // size of the QM zone
  this->natoms = qm_zone.qm.size() + qm_zone.link.size();

  // resize vectors
  this->attyp.resize(natoms);
  this->coord.resize(3 * natoms);

  // initialize XTB environment
  this->env = xtb_newEnvironment();
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }

  // initialize XTB calculator
  this->calc = xtb_newCalculator();
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }

  // initialize XTB results
  this->res = xtb_newResults();
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }

  // initialize QM atom types
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it, ++i) {
    this->attyp[i] = it->atomic_number;
  }
  // initialize QM links
  for (std::set<QM_Link>::const_iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it, ++i) {
    this->attyp[i] = it->atomic_number;
  }

  // process QM coordinates
  this->process_input_coordinates(topo, conf, sim, qm_zone);

  // Fortran API call - pass pointers
  this->mol = xtb_newMolecule(this->env, &(this->natoms), attyp.data(), coord.data(), 
                              &(this->charge), &(this->uhf), NULL, NULL);
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }

  // parametrize the molecule with GFN1-xTB or GFN2-xTB
  if (this->param->hamiltonian == 1) {
    xtb_loadGFN1xTB(this->env, this->mol, this->calc, NULL);
    if (xtb_checkEnvironment(this->env)) {
        xtb_showEnvironment(this->env, NULL);
        return 1;
    }
  }
  else if (this->param->hamiltonian == 2) {
    xtb_loadGFN2xTB(this->env, this->mol, this->calc, NULL);
    if (xtb_checkEnvironment(this->env)) {
        xtb_showEnvironment(this->env, NULL);
        return 1;
    }
  }

  // verbosity level
  xtb_setVerbosity(this->env, this->param->verbosity);
  if (xtb_checkEnvironment(env)) {
      xtb_showEnvironment(env, NULL);
      return 1;
  }

  DEBUG(15, "Initialized " << this->name());
  return 0;
}

int interaction::XTB_Worker::process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) {

  // first process QM coordinates
  this->process_input_coordinates(topo, conf, sim, qm_zone);

  // Fortran API call - pass pointers
  xtb_updateMolecule(this->env, this->mol, coord.data(), NULL);
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }

  // process point charges
  this->process_input_pointcharges(topo, conf, sim, qm_zone);

  // set external charges
  xtb_setExternalCharges(this->env, this->calc, &(this->num_charges), this->numbers.data(),
                         this->charges.data(), this->point_charges.data());
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }

  return 0;
}

void interaction::XTB_Worker::process_input_coordinates(const topology::Topology& topo
                                                      , const configuration::Configuration& conf
                                                      , const simulation::Simulation& sim
                                                      , const interaction::QM_Zone& qm_zone) {
  // MM -> QM length unit is inverse of input value
  double len_to_qm = 1.0 / this->param->unit_factor_length;

  // transfer QM coordinates
  DEBUG(15, "Transfering QM coordinates to XTB");
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->coord[i]     = it->pos(0) * len_to_qm;
    this->coord[i + 1] = it->pos(1) * len_to_qm;
    this->coord[i + 2] = it->pos(2) * len_to_qm;
    i += 3; // x, y, z component
  }
  // transfer capping atoms
  DEBUG(15, "Transfering capping atoms coordinates to XTB");
  for (std::set<QM_Link>::const_iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " "
      << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->coord[i]     = it->pos(0) * len_to_qm;
    this->coord[i + 1] = it->pos(1) * len_to_qm;
    this->coord[i + 2] = it->pos(2) * len_to_qm;
    i += 3; 
  } 
}

void interaction::XTB_Worker::process_input_pointcharges(const topology::Topology& topo
                                                       , const configuration::Configuration& conf
                                                       , const simulation::Simulation& sim
                                                       , const interaction::QM_Zone& qm_zone) {
  // MM -> QM length unit is inverse of input value
  double len_to_qm = 1.0 / this->param->unit_factor_length;
  double cha_to_qm = 1.0 / this->param->unit_factor_charge;

  unsigned int i = 0;
  for (std::set<MM_Atom>::const_iterator 
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    coord[i] = it->pos(0);
    coord[i + 1] = it->pos(1);
    coord[i + 2] = it->pos(2);
    i += 3; // x, y, z component
  }
  // scale coordinates
  std::transform(coord.begin(), coord.end(), coord.begin(), [len_to_qm](double& e) { return e *= len_to_qm; });
}

int interaction::XTB_Worker::run_calculation() {
  // run singlepoint calculation
  xtb_singlepoint(this->env, this->mol, this->calc, this->res);
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }
  // release charges after successful run
  xtb_releaseExternalCharges(this->env, this->calc);
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }
  return 0;
}

void interaction::XTB_Worker::print_coordinates(const std::vector<double>& coord) {
  for (unsigned int i = 0; i < coord.size(); i += 3) {
    std::cout << std::setw(25) << coord[i]
              << std::setw(25) << coord[i + 1]
              << std::setw(25) << coord[i + 2]
              << std::endl;
  }
}

int interaction::XTB_Worker::process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) {
  int err;

  // parse energy
  err = this->parse_energy(qm_zone);
  if (err) return err;

  // parse QM gradients
  err = this->parse_qm_gradients(qm_zone);
  if (err) return err;

  // parse MM gradients or charges
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // also parse MM gradients
    unsigned int ncharges = this->get_num_charges(sim, qm_zone);
    err = this->parse_mm_gradients(qm_zone, ncharges);
    if (err) return err;
  }
  else if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
             && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    // also extract charges
    err = this->parse_charges(qm_zone);
    if (err) return err;
  }

  return 0;
}

int interaction::XTB_Worker::parse_energy(interaction::QM_Zone& qm_zone) const {
  double energy; // energy will live here
  xtb_getEnergy(this->env, this->res, &energy);
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }
  qm_zone.QM_energy() = energy * this->param->unit_factor_energy;
  return 0;
}

int interaction::XTB_Worker::parse_qm_gradients(interaction::QM_Zone& qm_zone) const {
  std::vector<double> gradients(this->natoms * 3, 0.0); // gradients will live here
  xtb_getGradient(this->env, this->res, gradients.data());
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }
  unsigned int i = 0;
  for (std::set<QM_Atom>::iterator
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    it->force(0) = -this->param->unit_factor_force * gradients[i];
    it->force(1) = -this->param->unit_factor_force * gradients[i + 1];
    it->force(2) = -this->param->unit_factor_force * gradients[i + 2];
    i += 3;
  }
  i = 0;
  for (std::set<QM_Link>::iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    it->force(0) = -this->param->unit_factor_force * gradients[i];
    it->force(1) = -this->param->unit_factor_force * gradients[i + 1];
    it->force(2) = -this->param->unit_factor_force * gradients[i + 2];
    i += 3;
  }
  return 0;
}

int interaction::XTB_Worker::parse_mm_gradients(interaction::QM_Zone& qm_zone, const unsigned int ncharges) const {
  std::vector<double> pc_gradients(ncharges * 3, 0.0); // gradients will live here
  xtb_getPCGradient(this->env, this->res, pc_gradients.data());
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }
  unsigned int i = 0;
  for (std::set<MM_Atom>::iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    it->force(0) = -this->param->unit_factor_force * pc_gradients[i];
    it->force(1) = -this->param->unit_factor_force * pc_gradients[i + 1];
    it->force(2) = -this->param->unit_factor_force * pc_gradients[i + 2];
    i += 3;
  }
  return 0;
}

int interaction::XTB_Worker::parse_charges(interaction::QM_Zone& qm_zone) const {
  int err = 0;
  return err;
}