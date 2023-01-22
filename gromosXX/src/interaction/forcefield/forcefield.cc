/**
 * @file forcefield.cc
 * contains the inline functions for
 * forcefield.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../util/prepare_virial.h"

#include "forcefield.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

interaction::Forcefield::~Forcefield()
{
  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    delete *it;
  }
}

int interaction::Forcefield
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  
  int i = 0;

  for(iterator it = begin(), to = end();
      it != to;
      ++it){

    DEBUG(8, "init " << (*it)->name);
    i += (*it)->init(topo, conf, sim, os, quiet);
  }

  return i;
}



interaction::Interaction * interaction::Forcefield
::interaction(std::string name)
{
  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    if ((*it)->name == name) return *it;
  }
  return NULL;
}

interaction::Interaction const * interaction::Forcefield
::interaction(std::string name)const
{
  for(const_iterator it = begin(), to = end();
      it != to;
      ++it){
    if ((*it)->name == name) return *it;
  }
  return NULL;
}

int interaction::Forcefield
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim)
{
  DEBUG(5, "forcefield: calculate interaction");

  //m_timer.start(sim);

  conf.current().force = 0.0;
  
  if (sim.param().force.force_groups) {
    for(unsigned int i = 0; i < conf.special().force_groups.size(); ++i) {
      for(unsigned int j = 0; j < conf.special().force_groups.size(); ++j) {
        conf.special().force_groups[i][j] = 0.0;
      }
    }
  }
  
  DEBUG(15, "zero energies");
  conf.current().energies.zero();

  DEBUG(15, "zero lambda energies");
  conf.current().perturbed_energy_derivatives.zero();

  DEBUG(15, "zero sasa and sasavolume");
  conf.current().sasa_tot = 0.0;
  conf.current().sasa_buriedvol_tot = 0.0;

  conf.current().virial_tensor = 0.0;

  // prepare for the virial
  util::prepare_virial(topo, conf, sim);
  //if (sim.param().eds.eds) {//Todo: why no if? bschroed
  
    const unsigned int numstates = conf.special().eds.force_endstates.size();// sim.param().eds.numstates;
    DEBUG(15, "number of eds states " << numstates);
    assert(conf.special().eds.force_endstates.size() == numstates);
    assert(conf.special().eds.virial_tensor_endstates.size() == numstates);
    for (unsigned int state = 0; state < numstates; state++) {
      conf.special().eds.force_endstates[state] = 0.0;
      conf.special().eds.virial_tensor_endstates[state] = 0.0;
    }
 // }
  // ORIOL_GAMD
  if (sim.param().gamd.gamd) {
    DEBUG(15, "GAMD forcefield init"); 
    const unsigned int numaccelgroups = sim.param().gamd.igroups;
    DEBUG(15, "number of GAMD groups " << numaccelgroups);
    assert(conf.special().gamd.total_force.size() == numaccelgroups);
    assert(conf.special().gamd.dihe_force.size() == numaccelgroups);
    assert(conf.special().gamd.virial_tensor.size() == numaccelgroups);
    assert(conf.special().gamd.virial_tensor_dihe.size() == numaccelgroups);
    for (unsigned int group = 0; group < numaccelgroups; group++) {
      conf.special().gamd.dihe_force[group] = 0.0;
      conf.special().gamd.virial_tensor_dihe[group] = 0.0;
      conf.special().gamd.total_force[group] = 0.0;
      conf.special().gamd.virial_tensor[group] = 0.0;
    }
 }
 DEBUG(15, "Finished GAMD forcefield init");
  for (iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(5, "interaction: " << (*it)->name);
    // !!! crash if error
    int error=(*it)->calculate_interactions(topo, conf, sim);
    if (error){
      return 1;
    }
  }
  //m_timer.stop();  
  return 0;
}
    

void interaction::Forcefield
::print_timing(std::ostream & os)
{
  // m_timer.print(os);

  for(iterator it = begin(), to = end();
      it != to;
      ++it){

    (*it)->print_timing(os);

  }
}

