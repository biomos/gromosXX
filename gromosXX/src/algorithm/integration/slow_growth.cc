/**
 * @file slow_growth.cc
 * slow growth implementation
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

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
  // update masses and individual lambda-values
  topo.update_for_lambda();
  return 0;
}
