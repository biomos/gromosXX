/**
 * @file ghost_worker.cc
 * The worker class for the Ghost software (does nothing)
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
#include "ghost_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::Ghost_Worker::Ghost_Worker() : QM_Worker("Ghost Worker"), param(nullptr) {}; 

interaction::Ghost_Worker::~Ghost_Worker() = default;

int interaction::Ghost_Worker::init(const topology::Topology& topo
                                , const configuration::Configuration& conf
                                , simulation::Simulation& sim
                                , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());

  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.ghost);
  QM_Worker::param = this->param;

  // parent class initialization (trajectory files)
  int err = QM_Worker::init(topo, conf, sim, qm_zone);
  if (err) return err;

  DEBUG(15, "Initialized " << this->name());
  return 0;
}

int interaction::Ghost_Worker::process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) {
  int err = 0;
  // does nothing
  return err;
}

int interaction::Ghost_Worker::run_calculation() {
  DEBUG(15, "Ghost calculation");
  // does nothing
  return 0;
}

int interaction::Ghost_Worker::process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) {
  int err = 0;
  // does nothing
  return err;
}