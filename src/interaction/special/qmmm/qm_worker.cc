/**
 * @file qm_worker.cc
 * implements the factory function for the QM_Worker class
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../util/timing.h"

#include "qm_storage.h"
#include "qm_worker.h"
#include "mndo_worker.h"
#include "turbomole_worker.h"


interaction::QM_Worker * interaction::QM_Worker::get_instance(const simulation::Simulation & sim) {
  switch(sim.param().qmmm.software) {
    case simulation::qmmm_software_mndo :
      return new MNDO_Worker;
    case simulation::qmmm_software_turbomole :
      return new Turbomole_Worker;
    default:
      io::messages.add("QM worker not implemented", "QM_Worker", io::message::critical);
      break;
  }
  return NULL;
}
