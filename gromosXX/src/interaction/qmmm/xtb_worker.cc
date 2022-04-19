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

  // MM -> QM is inverse of input value
  double len_to_qm = 1.0 / this->param->unit_factor_length;

  // flatten the Vector3 to a Vector1 (which can be accessed as C style array by XTB)
  i = 0;
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    this->coord[i] = it->pos(0);
    this->coord[i + 1] = it->pos(1);
    this->coord[i + 2] = it->pos(2);
    i += 3; // x, y, z component
  }
  // scale coordinates
  std::transform(this->coord.begin(), this->coord.end(), this->coord.begin(), [len_to_qm](double& e) { return e *= len_to_qm; });

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

int interaction::XTB_Worker::run_calculation() {
  std::cout << "System call" << std::endl;
  return 0;
}

int interaction::XTB_Worker::process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) {


  std::cout << "Process input" << std::endl;

  return 0;
}

int interaction::XTB_Worker::process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) {
  std::cout << "Read output" << std::endl;
  return 0;
}