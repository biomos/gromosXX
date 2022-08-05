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

  // Charge and spin
  this->charge = qm_zone.charge() * this->param->unit_factor_charge;
  this->uhf = (qm_zone.spin_mult() - 1) / 2; // spin multiplicity to # unpaired electrons

  // size of the QM zone
  this->natoms = qm_zone.qm.size() + qm_zone.link.size();

  // resize vectors ahead of time
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

  // initialize QM atoms and QM link atoms
  this->init_input_atom_types(qm_zone);

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

void interaction::XTB_Worker::init_input_atom_types(const interaction::QM_Zone& qm_zone) {
  // QM atoms
  unsigned int i = 0;
  DEBUG(15, "Initializing QM atom types");
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it, ++i) {
    DEBUG(15, it->index << " " << it->atomic_number);
    this->attyp[i] = it->atomic_number;
  }
  // QM link atoms (iterator i keeps running)
  DEBUG(15, "Initializing capping atom types");
  for (std::set<QM_Link>::const_iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it, ++i) {
    DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " "
      << it->atomic_number);
    this->attyp[i] = it->atomic_number;
  }
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

  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // process point charges
    this->process_input_pointcharges(topo, conf, sim, qm_zone);

    if (this->ncharges > 0) { // set external charges only if > 0 ; otherwise xTB crashes
      xtb_setExternalCharges(this->env, this->calc, &(this->ncharges), this->numbers.data(),
                             this->charges.data(), this->point_charges.data());
      if (xtb_checkEnvironment(this->env)) {
          xtb_showEnvironment(this->env, NULL);
          return 1;
      }
    }
  }

  return 0;
}

void interaction::XTB_Worker::process_input_coordinates(const topology::Topology& topo
                                                      , const configuration::Configuration& conf
                                                      , const simulation::Simulation& sim
                                                      , const interaction::QM_Zone& qm_zone) {
  // Gromos -> QM length unit is inverse of input value from QM/MM specification file
  const double len_to_qm = 1.0 / this->param->unit_factor_length;

  // transfer QM coordinates
  DEBUG(15, "Transfering QM coordinates to XTB");
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    math::vector_c2f(this->coord, it->pos, i, len_to_qm);
    ++i;
  }
  // transfer capping atoms (index i keeps running...)
  DEBUG(15, "Transfering capping atoms coordinates to XTB");
  for (std::set<QM_Link>::const_iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " "
      << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    math::vector_c2f(this->coord, it->pos, i, len_to_qm);
    ++i;
  } 
}

void interaction::XTB_Worker::process_input_pointcharges(const topology::Topology& topo
                                                       , const configuration::Configuration& conf
                                                       , const simulation::Simulation& sim
                                                       , const interaction::QM_Zone& qm_zone) {
  // MM -> QM length unit is inverse of input value
  const double cha_to_qm = 1.0 / this->param->unit_factor_charge;
  const double len_to_qm = 1.0 / this->param->unit_factor_length;

  this->ncharges = this->get_num_charges(sim, qm_zone); // this also checks for COS
  this->numbers.resize(this->ncharges);
  this->charges.resize(this->ncharges);
  this->point_charges.resize(this->ncharges * 3);

  DEBUG(15, "Transfering point charges coordinates to XTB");

  unsigned int i = 0; // iterate over atoms - keep track of the offset for Fortran arrays
  for (std::set<MM_Atom>::const_iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    // memory layout of point charge arrays: 
    // numbers (1d), charges (1d), point_charges (3d): one-dimensional arrays 
    // COS numbers, charges, cartesian coordinates are after MM numbers, charges, cartesian coordinates
    if (it->is_polarisable) {
      // MM atom minus COS
      DEBUG(15, it->index << " " << it->atomic_number << " " 
        << (it->charge - it->cos_charge) * cha_to_qm << " " << math::v2s(it->pos * len_to_qm));
      this->numbers[i] = it->atomic_number;
      this->charges[i] = (it->charge - it->cos_charge) * cha_to_qm;
      math::vector_c2f(this->point_charges, it->pos, i, len_to_qm);
      ++i;
      // COS
      DEBUG(15, it->index << " " << it->atomic_number << " " 
        << it->cos_charge * cha_to_qm << " " << math::v2s((it->pos + it->cosV) * len_to_qm));
      this->numbers[i] = it->atomic_number;
      this->charges[i] = it->cos_charge * cha_to_qm;
      math::vector_c2f(this->point_charges, it->cosV, i, len_to_qm);
      ++i;
    }
    else {
      DEBUG(15, it->index << " " << it->atomic_number << " " 
        << it->charge * cha_to_qm << " " << math::v2s(it->pos * len_to_qm));
      this->numbers[i] = it->atomic_number;
      this->charges[i] = it->charge * cha_to_qm;
      math::vector_c2f(this->point_charges, it->pos, i, len_to_qm);
      ++i;
    }
  }
}

int interaction::XTB_Worker::run_calculation() {
  DEBUG(15, "Call XTB singlepoint calculation");
  // run singlepoint calculation
  xtb_singlepoint(this->env, this->mol, this->calc, this->res);
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }
  return 0;
}

int interaction::XTB_Worker::process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) {

  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical && this->ncharges > 0) {
    // release charges after successful run
    xtb_releaseExternalCharges(this->env, this->calc);
    if (xtb_checkEnvironment(this->env)) {
        xtb_showEnvironment(this->env, NULL);
        return 1;
    }
  }

  // parse energy
  int err = this->parse_energy(qm_zone);
  if (err) return err;

  // parse QM gradients
  err = this->parse_qm_gradients(qm_zone);
  if (err) return err;

  // parse MM gradients or charges
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical && this->ncharges > 0) {
    // also parse MM gradients
    err = this->parse_mm_gradients(qm_zone);
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
  const double force_to_mm = this->param->unit_factor_force;
  unsigned int i = 0;
  // Parse QM atoms
  for (std::set<QM_Atom>::iterator
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    // forces = negative gradient (!)
    DEBUG(15, "Parsing gradients of QM atom " << it->index);
    math::vector_f2c(gradients, it->force, i, -1.0 * force_to_mm);
    DEBUG(15, "Force: " << math::v2s(it->force));
    ++i;
  }
  // Parse capping atoms (index i keeps running...)
  for (std::set<QM_Link>::iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15, "Parsing gradient of capping atom " << it->qm_index << "-" << it->mm_index);
    math::vector_f2c(gradients, it->force, i, -1.0 * force_to_mm);
    DEBUG(15, "Force: " << math::v2s(it->force));
    ++i;
  }
  return 0;
}

int interaction::XTB_Worker::parse_mm_gradients(interaction::QM_Zone& qm_zone) const {
  std::vector<double> pc_gradients(this->ncharges * 3, 0.0); // gradients will live here
  xtb_getPCGradient(this->env, this->res, pc_gradients.data());
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }
  const double force_to_mm = this->param->unit_factor_force;
  // Parse MM atoms
  unsigned int i = 0;
  for (std::set<MM_Atom>::iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    // forces = negative gradient (!)
    DEBUG(15,"Parsing gradient of MM atom " << it->index);
    math::vector_f2c(pc_gradients, it->force, i, -1.0 * force_to_mm);
    DEBUG(15, "Force: " << math::v2s(it->force));
    if (it->is_polarisable) {
      ++i; // COS gradients live directly past the corresponding MM gradients
      DEBUG(15, "Parsing gradient of COS of MM atom " << it->index);
      math::vector_f2c(pc_gradients, it->cos_force, i, -1.0 * force_to_mm);
      DEBUG(15, "Force " << math::v2s(it->cos_force));
    }
    ++i;
  }
  return 0;
}

int interaction::XTB_Worker::parse_charges(interaction::QM_Zone& qm_zone) const {
  std::vector<double> charges(this->natoms, 0.0); // charges will live here
  xtb_getCharges(this->env, this->res, charges.data());
  if (xtb_checkEnvironment(this->env)) {
    xtb_showEnvironment(this->env, NULL);
    return 1;
  }
  const double cha_to_mm = this->param->unit_factor_charge;
  // QM atoms
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator
             it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; it++) {
    DEBUG(15, "Parsing charge of QM atom " << it->index);
    it->qm_charge = charges[i] * cha_to_mm;
    DEBUG(15, "Charge: " << it->qm_charge);
    ++i;
  }
  // Capping atoms (iterator i keeps running)
  for (std::set<QM_Link>::const_iterator
             it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Parsing charge of capping atom " << it->qm_index << "-" << it->mm_index);
    it->qm_charge = charges[i] * cha_to_mm;
    DEBUG(15, "Charge: " << it->qm_charge);
    ++i;
  }
  return 0;
}