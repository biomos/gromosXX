/**
 * @file forcefield.cc
 * contains the inline functions for
 * forcefield.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <util/prepare_virial.h>

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

  const double start = util::now();

  conf.current().force = 0.0;
  
  DEBUG(15, "zero energies");
  conf.current().energies.zero();

  DEBUG(15, "zero lambda energies");
  conf.current().perturbed_energy_derivatives.zero();

  conf.current().virial_tensor = 0.0;

  // prepare for the virial
  util::prepare_virial(topo, conf, sim);
  
  const unsigned int numstates = sim.param().eds.numstates;
  for(unsigned int state = 0; state < numstates; state++){
    conf.special().eds.force_endstates[state] = 0.0;
    conf.special().eds.virial_tensor_endstates[state] = 0.0;
  }
  
  for(iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    if ((*it)->calculate_interactions(topo, conf, sim))
      return 1;
  }
  
  if(sim.param().eds.eds){
    // interactions have been calculated - now apply eds Hamiltonian
    std::vector<long double> prefactors(numstates);
    
    // get beta
    assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
    const double beta = 1.0/(sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
    const double s = sim.param().eds.s;
    long double sum_prefactors = 0.0;
    const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
    const std::vector<double> & eir = sim.param().eds.eir;
    
    for(unsigned int state = 0; state < numstates; state++){
      const double pre = exp(-beta * s * (eds_vi[state] - eir[state]));
      prefactors[state] = pre;
      sum_prefactors += pre;
      DEBUG(7,"eds_vi[ " << state << "] = " << eds_vi[state]);
      DEBUG(7,"eir[" << state << "]" << eir[state]);
      DEBUG(7,"pre = " << pre);
    }
    // calculate eds Hamiltonian
    conf.current().energies.eds_vr = -1.0 / (beta * s) * log(sum_prefactors);
    DEBUG(7,"eds_vr = " << conf.current().energies.eds_vr);
    // calculate eds contribution ...
    const double sum_prefactors_i = 1.0 / sum_prefactors;
    
    for(unsigned int state = 0; state < numstates; state++){
      const long double pi = prefactors[state] * sum_prefactors_i;
      //std::cerr << "state = " << state << ", pi = " << pi << std::endl;
      // ... to forces
      for(unsigned int i = 0; i < topo.num_atoms(); i++ ) {
         conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
      }
      // ... to virial
      for (int a=0; a<3; ++a){
        for(int b=0; b<3; ++b){
          conf.current().virial_tensor(b, a) +=
          pi * conf.special().eds.virial_tensor_endstates[state](b, a);   
        }
      }
    } // loop over states
     
  }
 
  m_timing += util::now() - start;
  
  return 0;
}


void interaction::Forcefield
::print_timing(std::ostream & os)
{
  os << "    " 
     << std::setw(40) << std::left << name
     << std::setw(20) << m_timing << "\n";

  for(iterator it = begin(), to = end();
      it != to;
      ++it){

    (*it)->print_timing(os);

  }
}

