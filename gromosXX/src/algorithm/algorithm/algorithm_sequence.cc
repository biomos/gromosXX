/**
 * @file algorithm_sequence.cc
 * contains the inline functions for
 * algorithm_sequence.
 */

#include <util/stdheader.h>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>

#include "algorithm_sequence.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

algorithm::Algorithm_Sequence::Algorithm_Sequence()
  : std::vector<Algorithm *>()
{
}

algorithm::Algorithm_Sequence::~Algorithm_Sequence()
{
  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

int algorithm::Algorithm_Sequence
::run(topology::Topology & topo, 
      configuration::Configuration &conf,
      simulation::Simulation &sim)
{
  DEBUG(5, "Algorithm_Sequence: apply algorithm");

  int ret;
  
  for(Algorithm_Sequence::iterator 
	it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "algorithm: " << (*it)->name);
    if((ret = (*it)->apply(topo, conf, sim))){
      DEBUG(1, "ERROR in algorithm_sequence: bailing out!");
      return ret;
    }
  }
  return 0;
}