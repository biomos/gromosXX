
/**
 * @file eds.cc
 * contains the implementation
 * for the EDS class
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "eds.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * EDS step.
 */

int algorithm::EDS
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim)
 {
  m_timer.start();

  const unsigned int numstates = sim.param().eds.numstates;
  switch (sim.param().eds.form) {
    case simulation::aeds:
    case simulation::aeds_search_eir:
    case simulation::aeds_search_emax_emin:
    case simulation::aeds_search_all:
    {
      // interactions have been calculated - now apply eds Hamiltonian
      std::vector<double> prefactors(numstates);
      // get beta
      double beta = 0.0;
      if(!sim.param().stochastic.sd){
            assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
            beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
      }
      else
            beta = 1.0 / (sim.param().stochastic.temp * math::k_Boltzmann);
          
      std::vector<double> eds_vi = conf.current().energies.eds_vi;
      DEBUG(7, "eds_vi[0] = " << eds_vi[0]);
      DEBUG(7, "eds_vi[1] = " << eds_vi[1]);
      
      const std::vector<double> & eir = sim.param().eds.eir;
      unsigned int state_i = 0;
      unsigned int state_j = 1;
      double partA = -beta * (eds_vi[state_i] - eir[state_i]);
      double partB = -beta * (eds_vi[state_j] - eir[state_j]);
      DEBUG(7, "partA " << partA);
      DEBUG(7, "partB " << partB);
      double sum_prefactors = std::max(partA, partB)
        + log(1 + exp(std::min(partA, partB) - std::max(partA, partB)));
      prefactors[state_i] = partA;
      prefactors[state_j] = partB;

      for (unsigned int state = 2; state < numstates; state++) {
        double part = -beta * (eds_vi[state] - eir[state]);
        sum_prefactors = std::max(sum_prefactors, part)
          + log(1 + exp(std::min(sum_prefactors, part) - std::max(sum_prefactors, part)));
        prefactors[state] = part;
        DEBUG(7, "eds_vi[ " << state << "] = " << eds_vi[state]);
        DEBUG(7, "eir[" << state << "]" << eir[state]);
        DEBUG(7, "part = " << part);
      }

      // calculate eds Hamiltonian
      double demix = 0.0, kfac = 0.0, fkfac = 1.0;
      conf.current().energies.eds_vmix = -1.0 / beta * sum_prefactors;
      DEBUG(7, "eds_vmix = " << conf.current().energies.eds_vmix);

      // initilize search if necessary
      if (sim.param().eds.initaedssearch == true) {
        if (sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all) {
          sim.param().eds.emax = conf.current().energies.eds_vmix;
          sim.param().eds.emin = conf.current().energies.eds_vmix;
          sim.param().eds.searchemax = conf.current().energies.eds_vmix;
        }
        if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_all) {
          for (unsigned int is = 0; is < numstates; is++) {
            sim.param().eds.lnexpde[is] = (sim.param().eds.eir[is] - sim.param().eds.eir[0]) * -1.0 * beta;
          }
        }
      }

      // accelerate eds Hamiltonian
      if (conf.current().energies.eds_vmix <= sim.param().eds.emin) {
        conf.current().energies.eds_vr = conf.current().energies.eds_vmix;
      }
      else if (conf.current().energies.eds_vmix >= sim.param().eds.emax) {
        conf.current().energies.eds_vr = conf.current().energies.eds_vmix - 0.5 * (sim.param().eds.emax - sim.param().eds.emin);
      }
      else {
        demix = conf.current().energies.eds_vmix - sim.param().eds.emin;
        kfac = 1.0 / (sim.param().eds.emax - sim.param().eds.emin);
        fkfac = 1.0 - kfac * demix;
        conf.current().energies.eds_vr = conf.current().energies.eds_vmix - 0.5 * kfac * demix * demix;
      }

      // calculate eds contribution ...
      for (unsigned int state = 0; state < numstates; state++) {
        const long double pi = exp(prefactors[state] - sum_prefactors);
        //std::cerr << "state = " << state << ", pi = " << pi << std::endl;
        // ... to forces
        DEBUG(7, "prefactor = " << pi);
        for (unsigned int i = 0; i < topo.num_atoms(); i++) {
          if (conf.current().energies.eds_vmix <= sim.param().eds.emin || conf.current().energies.eds_vmix >= sim.param().eds.emax) {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i);
          }
          else {
            conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i) * fkfac;
          }
          DEBUG(9, "force current: " << i << " = " << math::v2s(conf.current().force(i)));
        }
        // ... to virial
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            if (conf.current().energies.eds_vmix <= sim.param().eds.emin || conf.current().energies.eds_vmix >= sim.param().eds.emax) {
              conf.current().virial_tensor(b, a) +=
                pi * conf.special().eds.virial_tensor_endstates[state](b, a);
            }
            else {
              conf.current().virial_tensor(b, a) +=
                pi * conf.special().eds.virial_tensor_endstates[state](b, a) * fkfac;
            }
          }
        }
      }

      // parameter search
      if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all) {
        DEBUG(7, "entering parameter search");
        // OFFSET search
        if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_all) {
          double tau = double(sim.param().eds.asteps) + double(sim.param().eds.bsteps - sim.param().eds.asteps) * double(sim.steps()) / double(sim.param().step.number_of_steps);
          double expde = 0.0, eiremin = 0.0, eiremax = 0.0, eirestar = 0.0, eirdemix = 0.0, eirkfac = 0.0;
          for (unsigned int is = 0; is < numstates; is++) {
            eiremin = sim.param().eds.emin + sim.param().eds.eir[is];
            eiremax = sim.param().eds.emax + sim.param().eds.eir[is];
            if (eds_vi[is] <= eiremin) {
              eirestar = eds_vi[is];
            }
            else if (eds_vi[is] >= eiremax) {
              eirestar = eds_vi[is] - 0.5 * (eiremax - eiremin);
            }
            else {
              eirdemix = eds_vi[is] - eiremin;
              eirkfac = 1.0 / (eiremax - eiremin);
              eirestar = eds_vi[is] - 0.5 * eirkfac * eirdemix * eirdemix;
            }
            expde = -1.0 * beta * (eirestar - conf.current().energies.eds_vr);
            sim.param().eds.lnexpde[is] += log(exp(-1.0 / tau) / (1.0 - exp(-1.0 / tau)));
            sim.param().eds.lnexpde[is] = std::max(sim.param().eds.lnexpde[is], expde) + log(1.0 + exp(std::min(sim.param().eds.lnexpde[is], expde) - std::max(sim.param().eds.lnexpde[is], expde)));
            sim.param().eds.lnexpde[is] += log(1.0 - exp(-1.0 / tau));
            sim.param().eds.statefren[is] = -1.0 / beta * sim.param().eds.lnexpde[is];
            sim.param().eds.eir[is] = sim.param().eds.statefren[is] - sim.param().eds.statefren[0];
          }
        }

        // EMAX and EMIN search
        if (sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all) {
          // find the state we are in
          unsigned int state = 0;
          double min = eds_vi[0] - eir[0];
          for (unsigned int is = 1; is < numstates; is++) {
            if ((eds_vi[is] - eir[is]) < min) {
              min = eds_vi[is] - eir[is];
              state = is;
            }
          }
          // mark the state we are in as visited
          sim.param().eds.visitedstates[state] = true;
          // monitor energy landscape
          sim.param().eds.visitcounts[state]++;
          double tempenergy = sim.param().eds.avgenergy[state] + (eds_vi[state] - sim.param().eds.avgenergy[state]) / double(sim.param().eds.visitcounts[state]);
          sim.param().eds.eiravgenergy[state] += (eds_vi[state] - sim.param().eds.eiravgenergy[state] - eir[state]) / double(sim.param().eds.visitcounts[state]);
          if (sim.param().eds.visitcounts[state] > 1) {
            sim.param().eds.bigs[state] += (eds_vi[state] - tempenergy) * (eds_vi[state] - sim.param().eds.avgenergy[state]);
            sim.param().eds.stdevenergy[state] = sqrt(sim.param().eds.bigs[state] / (double(sim.param().eds.visitcounts[state] - 1)));
          }
          sim.param().eds.avgenergy[state] = tempenergy;
          tempenergy = sim.param().eds.eiravgenergy[0];
          // find the state with the minimum average energy
          int targetstate = 0;
          for (unsigned int is = 1; is < numstates; is++) {
            if (tempenergy >= sim.param().eds.eiravgenergy[is] && sim.param().eds.visitcounts[is] > 1) {
              targetstate = is;
              tempenergy = sim.param().eds.eiravgenergy[is];
            }
          }
          double globminavg = sim.param().eds.eiravgenergy[targetstate];
          double globminfluc = sim.param().eds.stdevenergy[targetstate];
          // prevent rounding error in the first simulation step
          if (sim.param().eds.initaedssearch == true) {
            globminavg = sim.param().eds.emin;
          }
          // EMAX
          // search for the highest transition energy between states
          if (state != sim.param().eds.oldstate && sim.param().eds.searchemax < conf.current().energies.eds_vmix) {
            sim.param().eds.searchemax = conf.current().energies.eds_vmix;
          }
          // search for maximum transition energy for Emax in the beginning
          if (sim.param().eds.emaxcounts == 0) {
            sim.param().eds.emax = sim.param().eds.searchemax;
          }
          // check if we visited all the states; if so, stop updating searchemax, add it to Emax and reset
          unsigned int statecount = 0;
          for (unsigned int is = 0; is < numstates; is++) {
            if (sim.param().eds.visitedstates[is] == true) {
              statecount++;
            }
          }
          if (statecount == sim.param().eds.visitedstates.size()) {
            if (sim.param().eds.emaxcounts == 0) {
              sim.param().eds.emax = 0.0;
              sim.param().eds.fullemin = true;
            }
            sim.param().eds.emaxcounts += 1;
            sim.param().eds.emax += (sim.param().eds.searchemax - sim.param().eds.emax) / double(sim.param().eds.emaxcounts);
            sim.param().eds.searchemax = globminavg;
            for (unsigned int is = 0; is < numstates; is++) {
              sim.param().eds.visitedstates[is] = false;
            }
          }

          // EMIN
          double bmax = 0.0;
          if (sim.param().eds.bmaxtype == 1)
          {
            bmax = sim.param().eds.setbmax;
          }
          if (sim.param().eds.bmaxtype == 2)
          {
            // this implies that the fluctuations of the energies do not change with the acceleration
            bmax = sim.param().eds.stdevenergy[targetstate] * sim.param().eds.setbmax;
          }
          if ((sim.param().eds.emax - globminavg) <= bmax) {
            sim.param().eds.emin = sim.param().eds.emax;
            DEBUG(7, "emin1 " << sim.param().eds.emin);
          }
          else {
            sim.param().eds.emin = 2.0 * (globminavg + bmax) - sim.param().eds.emax;
            DEBUG(7, "emin2 " << sim.param().eds.emin);
            DEBUG(7, "emin2 " << globminavg);
            if (sim.param().eds.emin < globminavg) {
              sim.param().eds.emin = (-sim.param().eds.emax * sim.param().eds.emax
                + 2.0 * sim.param().eds.emax * bmax
                + 2.0 * sim.param().eds.emax * globminavg
                - globminavg * globminavg)
                / (2.0 * bmax);
              DEBUG(7, "emin3 " << sim.param().eds.emin);
            }
            // security measure to prevent extreme emins in the beginning of the simulation before we saw a full round-trip
            if (sim.param().eds.fullemin == false && sim.param().eds.emin < globminavg) {
              sim.param().eds.emin = globminavg;
              // globminavg can be larger than the current largest transition energy
              if (sim.param().eds.emin > sim.param().eds.emax) {
                sim.param().eds.emin = sim.param().eds.emax;
              }
              DEBUG(7, "emin4 " << sim.param().eds.emin);
            }
          }
          conf.current().energies.eds_globmin = globminavg;
          conf.current().energies.eds_globminfluc = globminfluc;
          sim.param().eds.oldstate = state;
        }
      }

      DEBUG(7, "updating energy configuration");

      conf.current().energies.eds_emax = sim.param().eds.emax;
      conf.current().energies.eds_emin = sim.param().eds.emin;
      for (unsigned int is = 0; is < numstates; is++) {
        conf.current().energies.eds_eir[is] = sim.param().eds.eir[is];
      }

      if (sim.param().eds.initaedssearch == true) {
        sim.param().eds.initaedssearch = false;
      }

      break;
    }
    case simulation::single_s:
    {
      // interactions have been calculated - now apply eds Hamiltonian
      std::vector<double> prefactors(numstates);
      // get beta
      double beta = 0.0;
      if(!sim.param().stochastic.sd){
            assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
            beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
      }
      else
            beta = 1.0 / (sim.param().stochastic.temp * math::k_Boltzmann);
          
      assert(sim.param().eds.s.size() == 1);
      const double s = sim.param().eds.s[0];
      
      std::vector<double> eds_vi = conf.current().energies.eds_vi;
      DEBUG(7, "eds_vi[0] = " << eds_vi[0]);
      DEBUG(7, "eds_vi[1] = " << eds_vi[1]);
     // DEBUG(7, "conf2 " << conf2);
      
      if (conf2 != NULL) {
        for (unsigned int i = 0; i < eds_vi.size(); ++i) {
          DEBUG(7, "conf2->current().energies.eds_vi[i] = " << conf2->current().energies.eds_vi[i]);
          eds_vi[i] += conf2->current().energies.eds_vi[i];
        }
      }
      
      const std::vector<double> & eir = sim.param().eds.eir;
      unsigned int state_i = 0;
      unsigned int state_j = 1;
      double partA = -beta * s * (eds_vi[state_i] - eir[state_i]);
      double partB = -beta * s * (eds_vi[state_j] - eir[state_j]);
      DEBUG(7, "partA " << partA);
      DEBUG(7, "partB " << partB);
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
    case simulation::multi_s:
    {
      // interactions have been calculated - now apply eds Hamiltonian
      std::vector<double> prefactors(numstates, 0.0);
      // get beta
      double beta = 0.0;
      if(!sim.param().stochastic.sd){
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
      }
      else
        beta = 1.0 / (sim.param().stochastic.temp * math::k_Boltzmann);
 
      const unsigned int numpairs = (numstates * (numstates - 1)) / 2;
      DEBUG(7, "number of eds states = " << numstates);
      DEBUG(7, "number of eds pairs = " << numpairs);
      assert(sim.param().eds.s.size() == numpairs);
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
              * (1 / s[pair_index] - 1);
      double pre_i_force = partA + elem2;
      double pre_j_force = partB + elem2;
      prefactors[state_i] = std::max(prefactors[state_i], pre_i_force)
              + log(1 + exp(std::min(prefactors[state_i], pre_i_force)));
      prefactors[state_j] = std::max(prefactors[state_j], pre_j_force)
              + log(1 + exp(std::min(prefactors[state_j], pre_j_force)
              - std::max(prefactors[state_j], pre_j_force)));
      DEBUG(7, "pre_i_force = " << pre_i_force << ", pre_j_force " << pre_j_force);
      DEBUG(7, "s_" << state_i << state_j << " = " << s[pair_index]);
      for (state_i = 1; state_i < (numstates - 1); state_i++) {
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
    case simulation::pair_s:
    {
      // interactions have been calculated - now apply eds Hamiltonian
      std::vector<double> prefactors(numstates, 0.0);
      // get beta
      double beta = 0.0;
      if(!sim.param().stochastic.sd){
        assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
        beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);
      }
      else
        beta = 1.0 / (sim.param().stochastic.temp * math::k_Boltzmann);

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
              * (1 / s[0] - 1);
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
  m_timer.stop();
  return 0;
} // eds

// m_timer.stop();

//  return 0;
//}


