/**
 * @file slow_growth.cc
 * slow growth implementation
 */

#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm.h>
#include "slow_growth.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * Slow growth update.
 */
int algorithm::Slow_Growth
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  // change the lambda value
  topo.lambda(topo.lambda() + sim.param().perturbation.dlamt
	      * sim.time_step_size());
  // update masses
  topo.update_for_lambda();
  return 0;
}
