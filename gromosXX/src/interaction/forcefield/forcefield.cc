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

  m_timer.start();

  conf.current().force = 0.0;
  
  DEBUG(15, "zero energies");
  conf.current().energies.zero();

  DEBUG(15, "zero lambda energies");
  conf.current().perturbed_energy_derivatives.zero();

  conf.current().virial_tensor = 0.0;

  // prepare for the virial
  util::prepare_virial(topo, conf, sim);
  //if (sim.param().eds.eds) {
  
    const unsigned int numstates = conf.special().eds.force_endstates.size();// sim.param().eds.numstates;
    DEBUG(15, "number of eds states " << numstates);
    assert(conf.special().eds.force_endstates.size() == numstates);
    assert(conf.special().eds.virial_tensor_endstates.size() == numstates);
    for (unsigned int state = 0; state < numstates; state++) {
      conf.special().eds.force_endstates[state] = 0.0;
      conf.special().eds.virial_tensor_endstates[state] = 0.0;
    }
 // }

  for (iterator it = begin(), to = end();
      it != to;
      ++it){
    DEBUG(7, "interaction: " << (*it)->name);
    if ((*it)->calculate_interactions(topo, conf, sim))
      return 1;
  }

  if (sim.param().eds.eds) {
    const unsigned int numstates = sim.param().eds.numstates;
    switch (sim.param().eds.form) {
      case simulation::single_s : {
        // interactions have been calculated - now apply eds Hamiltonian
        std::vector<long double> prefactors(numstates);
        // get beta
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
        assert(sim.param().eds.s.size()==1);
        const double s = sim.param().eds.s[0];
        long double sum_prefactors = 0.0;
        const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
        const std::vector<double> & eir = sim.param().eds.eir;

        for (unsigned int state = 0; state < numstates; state++) {
          const double pre = exp(-beta * s * (eds_vi[state] - eir[state]));
          prefactors[state] = pre;
          sum_prefactors += pre;
          DEBUG(7, "eds_vi[ " << state << "] = " << eds_vi[state]);
          DEBUG(7, "eir[" << state << "]" << eir[state]);
          DEBUG(7, "pre = " << pre);
        }
        // calculate eds Hamiltonian
        conf.current().energies.eds_vr = -1.0 / (beta * s) * log(sum_prefactors);
        DEBUG(7, "eds_vr = " << conf.current().energies.eds_vr);
        // calculate eds contribution ...
        const double sum_prefactors_i = 1.0 / sum_prefactors;

        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = prefactors[state] * sum_prefactors_i;
          //std::cerr << "state = " << state << ", pi = " << pi << std::endl;
          // ... to forces
          for (unsigned int i = 0; i < topo.num_atoms(); i++) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
          }
          // ... to virial
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              conf.current().virial_tensor(b, a) +=
                  pi * conf.special().eds.virial_tensor_endstates[state](b, a);
            }
          }
        } // loop over states

        break;
      }
      case simulation::multi_s : {
        // interactions have been calculated - now apply eds Hamiltonian
        std::vector<long double> prefactors(numstates,0.0);
        // get beta
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
        const unsigned int numpairs = (numstates * (numstates - 1)) / 2;
        DEBUG(7, "number of eds states = " << numstates);
        DEBUG(7, "number of eds pairs = " << numpairs);
        assert(sim.param().eds.s.size()==numpairs);
        const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
        const std::vector<double> & eir = sim.param().eds.eir;
        const std::vector<double> & s = sim.param().eds.s;
        long double sum_prefactors = 0.0;
        for (unsigned int state_i = 0; state_i < (numstates - 1); state_i++){
          for (unsigned int state_j = state_i + 1; state_j < numstates; state_j++){
            const unsigned int pair_index = (state_i*(2 * numstates - state_i-1)) / 2 + state_j - state_i - 1;
            DEBUG(7, "index of pair " << state_i << " - " << state_j << " = " << pair_index);
            // get the correct s_ij value
            // Hamiltonian
            const double s_ij = s[pair_index];
            assert(s_ij > 0);
            const double pre_i = exp(-beta * s_ij * (eds_vi[state_i] - eir[state_i]));
            const double pre_j = exp(-beta * s_ij * (eds_vi[state_j] - eir[state_j]));
            const double pre_i_pre_j = pre_i + pre_j;
            const double pipj_pow1 = pow(pre_i_pre_j,1.0/s_ij);
            const double pipj_pow2 = pow(pre_i_pre_j,1.0/s_ij - 1);
            sum_prefactors+=pipj_pow1;
            // additional calculation for forces
            const double pre_i_force = pipj_pow2*pre_i;
            const double pre_j_force = pipj_pow2*pre_j;
            prefactors[state_i] += pre_i_force;
            prefactors[state_j] += pre_j_force;
            DEBUG(7, "eds prefactor pair " << state_i << " - " << state_j << ", pre_i = " << pre_i << ", pre_j " << pre_j);
            DEBUG(7, "eds prefactor pair " << state_i << " - " << state_j << ", Vr component = " << pipj_pow1);
            DEBUG(7, "s_" << state_i << state_j << " = " << s_ij );
           
          }
        }
        // calculate eds Hamiltonian
        conf.current().energies.eds_vr = -1.0 / beta * log(sum_prefactors/(numstates-1));
        DEBUG(7, "eds_vr = " << conf.current().energies.eds_vr);
        // calculate eds contribution ...
        const double sum_prefactors_i = 1.0 / sum_prefactors;

        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = prefactors[state] * sum_prefactors_i;
          //std::cerr << "state = " << state << ", pi = " << pi << std::endl;
          // ... to forces
          for (unsigned int i = 0; i < topo.num_atoms(); i++) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
          }
          // ... to virial
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              conf.current().virial_tensor(b, a) +=
                  pi * conf.special().eds.virial_tensor_endstates[state](b, a);
            }
          }
        } // loop over states

        break;
      }
      case simulation::pair_s : {
        // interactions have been calculated - now apply eds Hamiltonian
        std::vector<long double> prefactors(numstates, 0.0);
        // get beta
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
        DEBUG(7, "number of eds states = " << numstates);
        assert(sim.param().eds.s.size() == numstates - 1);
        const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
        const std::vector<double> & eir = sim.param().eds.eir;
        const std::vector<double> & s = sim.param().eds.s;
        long double sum_prefactors = 0.0;
        // loop over (N-1) EDS state pairs
        for (unsigned int pair = 0; pair < sim.param().eds.pairs.size(); pair++) {
          assert(pair < s.size());
          const double s_ij = s[pair];
          assert(s_ij > 0);
          unsigned int state_i = sim.param().eds.pairs[pair].i - 1;
          unsigned int state_j = sim.param().eds.pairs[pair].j - 1;
          assert(state_i < eir.size() && state_i < eds_vi.size());
          assert(state_j < eir.size() && state_j < eds_vi.size());
          const double pre_i = exp(-beta * s_ij * (eds_vi[state_i] - eir[state_i]));
          const double pre_j = exp(-beta * s_ij * (eds_vi[state_j] - eir[state_j]));
          const double pre_i_pre_j = pre_i + pre_j;
          const double pipj_pow1 = pow(pre_i_pre_j, 1.0 / s_ij);
          const double pipj_pow2 = pow(pre_i_pre_j, 1.0 / s_ij - 1);
          sum_prefactors += pipj_pow1;
          // additional calculation for forces
          const double pre_i_force = pipj_pow2*pre_i;
          const double pre_j_force = pipj_pow2*pre_j;
          prefactors[state_i] += pre_i_force;
          prefactors[state_j] += pre_j_force;
          DEBUG(7, "eds prefactor pair " << state_i << " - " << state_j << ", pre_i = " << pre_i << ", pre_j " << pre_j);
          DEBUG(7, "eds prefactor pair " << state_i << " - " << state_j << ", Vr component = " << pipj_pow1);
          DEBUG(7, "s_" << state_i << state_j << " = " << s_ij);
        }
        // calculate eds Hamiltonian
        conf.current().energies.eds_vr = -1.0 / beta * log(sum_prefactors * numstates / ((numstates - 1)*2));
        DEBUG(7, "eds_vr = " << conf.current().energies.eds_vr);
        // calculate eds contribution ...
        const double sum_prefactors_i = 1.0 / sum_prefactors;

        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = prefactors[state] * sum_prefactors_i;
          // ... to forces
          for (unsigned int i = 0; i < topo.num_atoms(); i++) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
          }
          // ... to virial
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              conf.current().virial_tensor(b, a) +=
                  pi * conf.special().eds.virial_tensor_endstates[state](b, a);
            }
          }
        } // loop over states

        break;
      }
      default:
        io::messages.add("Unknown functional form of eds Hamiltonian. Should be 1 (single s), 2 (multi s), or 3 (N-1 pairs)",
            "Forcefield", io::message::critical);
    } // eds functional form switch
  } // eds
 
  m_timer.stop();
  
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

