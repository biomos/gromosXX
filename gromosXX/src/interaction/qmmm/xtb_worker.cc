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

  this->charge = qm_zone.charge() * this->param->unit_factor_charge;
  this->uhf = (qm_zone.spin_mult() - 1) / 2; // spin multiplicity to # unpaired electrons

  // size of the QM zone
  this->natoms = qm_zone.qm.size();

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

  // flatten the Vector3 to a Vector1 (which can be accessed as C style array by XTB)
  this->initialize_qm_coordinates(qm_zone);

  // initialize the molecule object
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

void interaction::XTB_Worker::initialize_qm_coordinates(const interaction::QM_Zone& qm_zone) {
  // MM -> QM length unit is inverse of input value
  double len_to_qm = 1.0 / this->param->unit_factor_length;
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    this->coord[i] = it->pos(0);
    this->coord[i + 1] = it->pos(1);
    this->coord[i + 2] = it->pos(2);
    i += 3; // x, y, z component
  }
  // scale coordinates
  std::transform(this->coord.begin(), this->coord.end(), this->coord.begin(), [len_to_qm](double& e) { return e *= len_to_qm; });
}

void interaction::XTB_Worker::initialize_mm_coordinates(const interaction::QM_Zone& qm_zone
                                                      , std::vector<double>& coord) {
  // MM -> QM length unit is inverse of input value
  double len_to_qm = 1.0 / this->param->unit_factor_length;
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

int interaction::XTB_Worker::process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) {

  // number of point charges
  int num_charges = this->get_num_charges(sim, qm_zone);

  // charges of the point charges
  std::vector<double> charges(num_charges, 0.0);

  // definition of chemical hardness parameters
  // for gradients on point charges:
  // https://xtb-docs.readthedocs.io/en/latest/pcem.html
  std::vector<int> numbers(num_charges, 0);

  // coordinates of point charges
  std::vector<double> point_charges(num_charges * 3, 0.0);

  // initialize QM atom coordinates
  this->initialize_qm_coordinates(qm_zone);

  xtb_updateMolecule(this->env, this->mol, this->coord.data(), NULL);
  if (xtb_checkEnvironment(this->env)) {
      xtb_showEnvironment(this->env, NULL);
      return 1;
  }

  // initialize MM atom charges and chemical hardness
  unsigned int i = 0;
  for (std::set<MM_Atom>::const_iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it, ++i) {
    numbers[i] = 7; // it seems like that this is a chemical hardness ...
    charges[i] = it->charge;
  }

  // initialize MM atom coordintes
  this->initialize_mm_coordinates(qm_zone, point_charges);

  // set external charges
  xtb_setExternalCharges(this->env, this->calc, &num_charges, numbers.data(),
                         charges.data(), point_charges.data());
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