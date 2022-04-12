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

int interaction::XTB_Worker::init(simulation::Simulation& sim) {
  DEBUG(15, "Initializing " << this->name());

  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.xtb);
  QM_Worker::param = this->param;

  xtb_newEnvironment();

  DEBUG(15, "Initialized " << this->name());
  return 0;
}

int interaction::XTB_Worker::system_call() {
  std::cout << "System call" << std::endl;
  return 0;
}

int interaction::XTB_Worker::write_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) {
  std::cout << "Write input" << std::endl;
  return 0;
}

int interaction::XTB_Worker::read_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) {
  std::cout << "Read output" << std::endl;
  return 0;
}