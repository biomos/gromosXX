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

  m_timer.start();

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
    DEBUG(5, "interaction: " << (*it)->name);
    if ((*it)->calculate_interactions(topo, conf, sim))
      return 1;
    DEBUG(5, "force old = " << math::v2s(conf.old().force(2)));
    DEBUG(5, "force current = " << math::v2s(conf.current().force(2)));
    if (sim.param().eds.eds)
      DEBUG(5, "force special = " << math::v2s(conf.special().eds.force_endstates[0](2)));
  }

  if (sim.param().eds.eds) {
    const unsigned int numstates = sim.param().eds.numstates;
    switch (sim.param().eds.form) {
      case simulation::single_s : {
        // interactions have been calculated - now apply eds Hamiltonian
        std::vector<double> prefactors(numstates);
        // get beta
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
        assert(sim.param().eds.s.size()==1);
        const double s = sim.param().eds.s[0];
        const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
        const std::vector<double> & eir = sim.param().eds.eir;
        unsigned int state_i = 0;
        unsigned int state_j = 1;
        double partA = -beta * s * (eds_vi[state_i] - eir[state_i]);
        double partB = -beta * s * (eds_vi[state_j] - eir[state_j]);
        double sum_prefactors = std::max(partA, partB)
                + log(1 + exp(std::min(partA, partB) - std::max(partA, partB)));
        prefactors[state_i] = partA;
        prefactors[state_j] = partB;

        for (unsigned int state = 2; state < numstates; state++) {
          double part = -beta * s * (eds_vi[state] - eir[state]);
          sum_prefactors = std::max(sum_prefactors, part)
                    + log(1 + exp(std::min(sum_prefactors, part) - std::max(sum_prefactors, part)));
          prefactors[state] = part;
          DEBUG(7, "eds_vi[ " << state << "] = " << eds_vi[state]);
          DEBUG(7, "eir[" << state << "]" << eir[state]);
          DEBUG(7, "part = " << part);
        }
        // calculate eds Hamiltonian
        conf.current().energies.eds_vr = -1.0 / (beta * s) * sum_prefactors;
        DEBUG(7, "eds_vr = " << conf.current().energies.eds_vr);

        // calculate eds contribution ...
        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = exp(prefactors[state] - sum_prefactors);
          //std::cerr << "state = " << state << ", pi = " << pi << std::endl;
          // ... to forces
          DEBUG(7, "prefactor = " << pi);
          for (unsigned int i = 0; i < topo.num_atoms(); i++) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
            DEBUG(9, "force current: " << i << " = " << math::v2s(conf.current().force(i)));
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
        std::vector<double> prefactors(numstates,0.0);
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
        // first pair
        unsigned int state_i = 0;
        unsigned int state_j = 1;
        unsigned int pair_index = 0;
        double partA = -beta * s[pair_index] * (eds_vi[state_i] - eir[state_i]);
        double partB = -beta * s[pair_index] * (eds_vi[state_j] - eir[state_j]);
        double sum_prefactors = (std::max(partA, partB)
                + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                / s[pair_index];
        double elem2 = (std::max(partA, partB)
                  + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                  * (1/s[pair_index] -1);
        double pre_i_force = partA + elem2;
        double pre_j_force = partB + elem2;
        prefactors[state_i] = std::max(prefactors[state_i], pre_i_force)
                  + log(1 + exp(std::min(prefactors[state_i], pre_i_force)));
        prefactors[state_j] = std::max(prefactors[state_j], pre_j_force)
                  + log(1 + exp(std::min(prefactors[state_j], pre_j_force)
                  - std::max(prefactors[state_j], pre_j_force)));
        DEBUG(7, "pre_i_force = " << pre_i_force << ", pre_j_force " << pre_j_force);
        DEBUG(7, "s_" << state_i << state_j << " = " << s[pair_index]);
        for (state_i = 1; state_i < (numstates - 1); state_i++){
          for (state_j = state_i + 1; state_j < numstates; state_j++) {
            pair_index = (state_i * (2 * numstates - state_i - 1)) / 2 + state_j - state_i - 1;
            DEBUG(7, "index of pair " << state_i << " - " << state_j << " = " << pair_index);
            // get the correct s_ij value
            const double s_ij = s[pair_index];
            assert(s_ij > 0);
            partA = -beta * s_ij * (eds_vi[state_i] - eir[state_i]);
            partB = -beta * s_ij * (eds_vi[state_j] - eir[state_j]);
            double elem = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    / s_ij;
            elem2 = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    * (1 / s_ij - 1);
            sum_prefactors = std::max(sum_prefactors, elem)
                    + log(1 + exp(std::min(sum_prefactors, elem) - std::max(sum_prefactors, elem)));
            DEBUG(7, "partA = " << partA << " , partB = " << partB << " , elem = "
                    << elem << " , sum = " << sum_prefactors);
            // additional calculation for forces
            pre_i_force = partA + elem2;
            pre_j_force = partB + elem2;
            prefactors[state_i] = std::max(prefactors[state_i], pre_i_force)
                    + log(1 + exp(std::min(prefactors[state_i], pre_i_force)
                    - std::max(prefactors[state_i], pre_i_force)));
            prefactors[state_j] = std::max(prefactors[state_j], pre_j_force)
                    + log(1 + exp(std::min(prefactors[state_j], pre_j_force)
                    - std::max(prefactors[state_j], pre_j_force)));
            DEBUG(7, "pre_i_force = " << pre_i_force << ", pre_j_force " << pre_j_force);
            DEBUG(7, "s_" << state_i << state_j << " = " << s_ij);
          }
        }
        // calculate eds Hamiltonian
        conf.current().energies.eds_vr = -1.0 / beta * (sum_prefactors + log(1.0 / double(numstates - 1)));
        DEBUG(7, "eds_vr = " << conf.current().energies.eds_vr);

        // calculate eds contribution ...
        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = exp(prefactors[state] - sum_prefactors);
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
        std::vector<double> prefactors(numstates, 0.0);
        // get beta
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
        DEBUG(7, "number of eds states = " << numstates);
        assert(sim.param().eds.s.size() == numstates - 1);
        const std::vector<double> & eds_vi = conf.current().energies.eds_vi;
        const std::vector<double> & eir = sim.param().eds.eir;
        const std::vector<double> & s = sim.param().eds.s;
        unsigned int state_i = sim.param().eds.pairs[0].i - 1;
        unsigned int state_j = sim.param().eds.pairs[0].j - 1;
        double partA = -beta * s[0] * (eds_vi[state_i] - eir[state_i]);
        double partB = -beta * s[0] * (eds_vi[state_j] - eir[state_j]);
        double sum_prefactors = (std::max(partA, partB)
                + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                / s[0];
        double elem2 = (std::max(partA, partB)
                  + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                  * (1/s[0] -1);
        double pre_i_force = partA + elem2;
        double pre_j_force = partB + elem2;
        prefactors[state_i] = std::max(prefactors[state_i], pre_i_force)
                  + log(1 + exp(std::min(prefactors[state_i], pre_i_force)));
        prefactors[state_j] = std::max(prefactors[state_j], pre_j_force)
                  + log(1 + exp(std::min(prefactors[state_j], pre_j_force)
                  - std::max(prefactors[state_j], pre_j_force)));
        DEBUG(7, "pre_i_force = " << pre_i_force << ", pre_j_force " << pre_j_force);
        DEBUG(7, "s_" << state_i << state_j << " = " << s[0]);

        // loop over (N-1) EDS state pairs
        for (unsigned int pair = 1; pair < sim.param().eds.pairs.size(); pair++) {
          assert(pair < s.size());
          const double s_ij = s[pair];
          assert(s_ij > 0);
          unsigned int state_i = sim.param().eds.pairs[pair].i - 1;
          unsigned int state_j = sim.param().eds.pairs[pair].j - 1;
          assert(state_i < eir.size() && state_i < eds_vi.size());
          assert(state_j < eir.size() && state_j < eds_vi.size());
          partA = -beta * s_ij * (eds_vi[state_i] - eir[state_i]);
          partB = -beta * s_ij * (eds_vi[state_j] - eir[state_j]);
          double elem = (std::max(partA, partB)
                  + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                  / s_ij;
          elem2 = (std::max(partA, partB)
                  + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                  * (1/s_ij -1);
          sum_prefactors = std::max(sum_prefactors, elem)
                  + log(1 + exp(std::min(sum_prefactors, elem) - std::max(sum_prefactors, elem)));
          DEBUG(7, "partA = " << partA << " , partB = " << partB << " , elem = "
                  << elem << " , sum = " << sum_prefactors);
          // additional calculation for forces
          pre_i_force = partA + elem2;
          pre_j_force = partB + elem2;
          prefactors[state_i] = std::max(prefactors[state_i], pre_i_force)
                  + log(1 + exp(std::min(prefactors[state_i], pre_i_force)
                  - std::max(prefactors[state_i], pre_i_force)));
          prefactors[state_j] = std::max(prefactors[state_j], pre_j_force)
                  + log(1 + exp(std::min(prefactors[state_j], pre_j_force)
                  - std::max(prefactors[state_j], pre_j_force)));
          DEBUG(7, "pre_i_force = " << pre_i_force << ", pre_j_force " << pre_j_force);
          DEBUG(7, "s_" << state_i << state_j << " = " << s_ij);
        }
        // calculate eds Hamiltonian
        DEBUG(7, "sum = " << sum_prefactors << " , num_states = " << numstates);
        conf.current().energies.eds_vr = -1.0 / beta * (sum_prefactors + log(double(numstates) / double((numstates - 1)*2)));
        DEBUG(3, "eds_vr = " << conf.current().energies.eds_vr);

        // calculate eds contribution ...
        for (unsigned int state = 0; state < numstates; state++) {
          const long double pi = exp(prefactors[state] - sum_prefactors);
          DEBUG(3, "pi = " << pi << " , prefactors = " << prefactors[state]);
          // ... to forces
          for (unsigned int i = 0; i < topo.num_atoms(); i++) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
            DEBUG(7, "force = " << math::v2s(conf.special().eds.force_endstates[state](i)));
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

