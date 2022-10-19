
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

#include "../../util/error.h"

#include "eds.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration




// gamma acceleration functions
// private functions

/**finds the gamma power needed to fulfill the target_frac */
void algorithm::EDS::find_gamma_pow(){
    double temp_gamma_pow = (1. / target_frac - 1) / 2;
    gamma_pow =  2*ceil(temp_gamma_pow);
    gamma_pow_1 = gamma_pow + 1;
    gamma_c = 1. / (2*(gamma_pow_1));
}

/**
 * calculates the linear (and the gamma term) scaling factor to get the target_frac
 * sets lin_frac, gamma_frac and gamma_pow based on target_frac
 */
void algorithm::EDS::find_lin_scaling_frac(){
    find_gamma_pow();
    double gama_frac_min = 1. / (gamma_pow + 1);
    lin_frac = (target_frac - gama_frac_min) / (1. - gama_frac_min);
    gamma_frac = 1. - lin_frac;
    // new transformation function (derivative of it) is
    // gamma_frac * gamma_pow_term + (1-gamma_frac) * gamma_pow_term
}

/**
* calculates f(x) for a given gamma power
* f(x) i the scaling factor for the energy (acceleration)
* where x = (H - Emin) / (Emax - Emin)
* gamma_pow has to be 0, 2, 4, 6...
* @param x 
*/
void algorithm::EDS::f_k_pow(double x){
    f_gamma = pow(2.0, gamma_pow) * pow(x - 0.5, gamma_pow_1) / gamma_pow_1 + gamma_c;
}

/**
 * calculates f'(x) for a given gamma power
 * f'(x) i the scaling factor for the force (acceleration)
 * where x = (H - Emin) / (Emax - Emin)
 * gamma_pow has to be 0, 2, 4, 6...
 * @param x 
 */
void algorithm::EDS::f_der_k_pow(double x){
    f_der_gamma = pow(2 * x - 1, gamma_pow);
}

/**
 * the 2 functions from above in one (optimized version)
 * @param x 
 */
void algorithm::EDS::f_f_der_k_pow(double x){
    f_der_gamma = pow(2 * x - 1, gamma_pow); // same as pow(2.0, gamma_pow) * pow(x - 0.5, gamma_pow)
    f_gamma = f_der_gamma * (x - 0.5) / gamma_pow_1 + gamma_c;
}

// calculate f(x) and f'(x) needed to get H* and F*

/** 
 * calculate f(x) needed to get H*
 * @param x 
 * @return f(x)
 */
double algorithm::EDS::get_f_x_lin(double x){
    f_k_pow(x);
    return gamma_frac * f_gamma + lin_frac * x;
}
/** 
 * calculate f'(x) needed to get F*
 * @param x 
 * @return f'(x)
 */
double algorithm::EDS::get_f_der_x_lin(double x){
    f_der_k_pow(x);
    return gamma_frac * f_der_gamma + lin_frac;
}

/**
 * calculate both f(x) and f'(x) at the same time
 * @param x 
 * @param f_x pointer to a double to assign the value of f(x)
 * @param f_der_x pointer to a double to assign the value of f'(x)
 */
void algorithm::EDS::get_f_f_der_x_lin(double x, double *f_x, double *f_der_x){
    f_f_der_k_pow(x);
    *f_x = gamma_frac * f_gamma + lin_frac * x;
    *f_der_x = gamma_frac * f_der_gamma + lin_frac;
}

// transition from f(x) to accelerated H

/**get the target_frac based on emin, emax and target_emax*/
void algorithm::EDS::get_target_frac(){
    if (emax == emin)
      target_frac = 1.1;
    else
      target_frac = (target_emax - emin) / (emax - emin);
}


// public functions

/**
 * sets emin, emax and target_emax 
 * also updates other gamma acceleration variables (e.g. target_frac, diff_emm, gamma_pow, ...)
 * @param emin_value
 * @param emax_value
 * @param target_emax_value
 */
void algorithm::EDS::set_EDS_params(double emin_value, double emax_value, double target_emax_value){
    emin = emin_value;
    emax = emax_value;
    target_emax = target_emax_value;
    update_gamma_params();
}

/** updates gamma acceleration variables (e.g. target_frac, diff_emm, gamma_pow, ...)*/
void algorithm::EDS::update_gamma_params(){
    diff_emm = emax - emin;
    get_target_frac();
    find_lin_scaling_frac();
}

/**
 * transforms energy value into x variable for gamma acceleration
 * @param E
 */
double algorithm::EDS::transform_E_2_x(double E){
    return (E - emin) / diff_emm;
}

/**
 * calculate accelerated energy & scaling factor for the force (gamma acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param f_der_x pointer to a variable to assign the value of the force scaling factor
 * @param E energy to be accelerated
 * @param E_offset energy offset (assumed 0 if not given)
 */
void algorithm::EDS::accelerate_E_F_gamma(double *E_a, double *f_der_x, double E, double E_offset=0){
    *f_der_x = 1.;
    if (target_emax >= emax)
        *E_a = E;
    double E_with_offset = E - E_offset;
    if (E_with_offset <= emin)
        *E_a = E;
    else if (E_with_offset >= emax)
        *E_a = target_emax + E - emax;
    else{
        double x = transform_E_2_x(E_with_offset);
        double f_x;
        get_f_f_der_x_lin(x, &f_x, f_der_x);
        *E_a = emin + f_x * (E_with_offset - emin) + E_offset;
    }
}

/**
 * calculate accelerated energy (gamma acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param E energy to be accelerated
 * @param E_offset energy offset (assumed 0 if not given)
 */
void algorithm::EDS::accelerate_E_gamma(double *E_a, double E, double E_offset=0){
    if (target_emax >= emax)
        *E_a = E;
    double E_with_offset = E - E_offset;
    if (E_with_offset <= emin)
        *E_a = E;
    else if (E_with_offset >= emax)
        *E_a = target_emax + E - emax;
    else{
        double x = transform_E_2_x(E_with_offset);
        double f_x = get_f_x_lin(x);
        *E_a = emin + f_x * (E_with_offset - emin) + E_offset;
    }
}


/**
 * calculate accelerated energy & scaling factor for the force (gaussian acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param fkfac pointer to a variable to assign the value of the force scaling factor
 * @param E energy to be accelerated
 * @param E_offset energy offset (assumed 0 if not given)
 */
void algorithm::EDS::accelerate_E_F_gauss(double *E_a, double *fkfac, double E, double E_offset=0){
    *fkfac = 1.;
    double E_with_offset = E - E_offset;
    if (E_with_offset <= emin)
        *E_a = E;
    else if (E_with_offset >= emax)
        *E_a = E - 0.5 * emax - emin;
    else{
        double demix = E_with_offset - emin;
        double kfac = 1.0 / diff_emm;
        *fkfac = 1.0 - kfac * demix;
        *E_a = E - 0.5 * kfac * demix * demix;
    }
}

/**
 * calculate accelerated energy & scaling factor for the force (gaussian acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param E energy to be accelerated
 * @param E_offset energy offset (assumed 0 if not given)
 */
void algorithm::EDS::accelerate_E_gauss(double *E_a, double E, double E_offset=0){
    double E_with_offset = E - E_offset;
    if (E_with_offset <= emin)
        *E_a = E;
    else if (E_with_offset >= emax)
        *E_a = E - 0.5 * emax - emin;
    else{
        double demix = E_with_offset - emin;
        double kfac = 1.0 / diff_emm;
        *E_a = E - 0.5 * kfac * demix * demix;
    }
}



/**
 * EDS step
 */

int algorithm::EDS
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim)
 {
  m_timer.start();

  int flag_acceleration_type = 1; // gamma acceleration
  // int flag_acceleration_type = 0; // old/gauss acceleration
  // this should be implemented through imd file...

  const unsigned int numstates = sim.param().eds.numstates;
  switch (sim.param().eds.form) {
    case simulation::aeds:
    case simulation::aeds_search_eir:
    case simulation::aeds_search_emax_emin:
    case simulation::aeds_search_all:
    case simulation::aeds_advanced_search:
    case simulation::aeds_advanced_search2:
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
        if (sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all  
        || sim.param().eds.form == simulation::aeds_advanced_search || sim.param().eds.form == simulation::aeds_advanced_search2) {
          sim.param().eds.emax = conf.current().energies.eds_vmix;
          sim.param().eds.emin = conf.current().energies.eds_vmix;
          sim.param().eds.searchemax = conf.current().energies.eds_vmix;
        }
        if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_all
         || sim.param().eds.form == simulation::aeds_advanced_search || sim.param().eds.form == simulation::aeds_advanced_search2) {
          for (unsigned int is = 0; is < numstates; is++) {
            sim.param().eds.lnexpde[is] = (sim.param().eds.eir[is] - sim.param().eds.eir[0]) * -1.0 * beta;
          }
        }
      }

      // update EDS params and call acceleration function 
      // update EDS params
      set_EDS_params(sim.param().eds.emin, sim.param().eds.emax, sim.param().eds.target_emax);
      // accelerate eds Hamiltonian
      if (flag_acceleration_type == 0)
        accelerate_E_F_gauss(&conf.current().energies.eds_vr, &fkfac, conf.current().energies.eds_vmix);
      else if (flag_acceleration_type == 1)
        accelerate_E_F_gamma(&conf.current().energies.eds_vr, &fkfac, conf.current().energies.eds_vmix);

      // calculate eds contribution ...
      for (unsigned int state = 0; state < numstates; state++) {
        const long double pi = exp(prefactors[state] - sum_prefactors);
        // ... to forces
        DEBUG(7, "prefactor = " << pi);
        for (unsigned int i = 0; i < topo.num_atoms(); i++) {
          // fkfac = 1. if outside of acceleration range
          conf.current().force(i) += pi * conf.special().eds.force_endstates[state](i) * fkfac;
          DEBUG(9, "force current: " << i << " = " << math::v2s(conf.current().force(i)));
        }
        // ... to virial
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            conf.current().virial_tensor(b, a) +=
              pi * conf.special().eds.virial_tensor_endstates[state](b, a) * fkfac;
          }
        }
      }

      // parameter search
      if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_emax_emin || 
          sim.param().eds.form == simulation::aeds_search_all || sim.param().eds.form == simulation::aeds_advanced_search
          || sim.param().eds.form == simulation::aeds_advanced_search2) {
        DEBUG(7, "entering parameter search");
        // OFFSET search
        if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_all 
        ||  sim.param().eds.form == simulation::aeds_advanced_search ||  sim.param().eds.form == simulation::aeds_advanced_search2) {
          double tau = double(sim.param().eds.asteps) + double(sim.param().eds.bsteps - sim.param().eds.asteps) * double(sim.steps()) / double(sim.param().step.number_of_steps);
          double expde = 0.0, eiremin = 0.0, eiremax = 0.0, eirestar = 0.0, eirdemix = 0.0, eirkfac = 0.0;
          for (unsigned int is = 0; is < numstates; is++) {
            // get the acceleration
            if (flag_acceleration_type == 0)
              accelerate_E_gauss(&eirestar, eds_vi[is]);
            else if (flag_acceleration_type == 1)
              accelerate_E_gamma(&eirestar, eds_vi[is]);
            
            // if in conventional search algorithm recalculate offsets using time decay function and update them
            if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_all){
              expde = -1.0 * beta * (eirestar - conf.current().energies.eds_vr);
              sim.param().eds.lnexpde[is] += log(exp(-1.0 / tau) / (1.0 - exp(-1.0 / tau)));
              sim.param().eds.lnexpde[is] = std::max(sim.param().eds.lnexpde[is], expde) + log(1.0 + exp(std::min(sim.param().eds.lnexpde[is], expde) - std::max(sim.param().eds.lnexpde[is], expde)));
              sim.param().eds.lnexpde[is] += log(1.0 - exp(-1.0 / tau));
              sim.param().eds.statefren[is] = -1.0 / beta * sim.param().eds.lnexpde[is];
              sim.param().eds.eir[is] = sim.param().eds.statefren[is] - sim.param().eds.statefren[0];
            }
            // if in adaptive fast search algorithm, calculate the free energies from the endstates to the reference without using time decay
            // and without updating the offsets
            if (sim.param().eds.form == simulation::aeds_advanced_search || sim.param().eds.form == simulation::aeds_advanced_search2){
              expde = -1.0 * beta * (eirestar - conf.current().energies.eds_vr);
              sim.param().eds.lnexpde[is] += log(exp(-1.0 / tau) / (1.0 - exp(-1.0 / tau)));
              sim.param().eds.lnexpde[is] = std::max(sim.param().eds.lnexpde[is], expde) + log(1.0 + exp(std::min(sim.param().eds.lnexpde[is], expde) - std::max(sim.param().eds.lnexpde[is], expde)));
              sim.param().eds.lnexpde[is] += log(1.0 - exp(-1.0 / tau));
              sim.param().eds.statefren[is] = -1.0 / beta * sim.param().eds.lnexpde[is];

            }
          }
        }

        // EMAX and EMIN search
        if (sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all || sim.param().eds.form == simulation::aeds_advanced_search) {
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

          
          sim.param().eds.target_emax = bmax;
          if (flag_acceleration_type == 0){
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
          }
          else if (flag_acceleration_type == 1){
            sim.param().eds.target_emax = bmax;
            sim.param().eds.emin = globminavg;
          }


          conf.current().energies.eds_globmin = globminavg;
          conf.current().energies.eds_globminfluc = globminfluc;
          sim.param().eds.oldstate = state;
        }

        // Adaptive AEDS search
        // at this point we have already computed the free energies of the difference endstates to the reference
        if (sim.param().eds.form == simulation::aeds_advanced_search || sim.param().eds.form == simulation::aeds_advanced_search2){
          // compute the prevalence of each endstate over this frame
          std::vector<double> exp_vi(numstates);
          double prevalence = 0.0;
          double total_endstates_ene = 0.0;
          for (unsigned int state = 0; state < numstates; state++) {
            exp_vi[state] = exp(-1.0 * beta * (eds_vi[state] - eir[state]));
            total_endstates_ene += exp_vi[state];
          }

          for (unsigned int state = 0; state < numstates; state++) {
            prevalence = exp_vi[state]/ total_endstates_ene;
            // Now check if the frame is contributing to the free energy
            // only states that have contributed at least 1/numstates on this frame will be checked to avoid counting frames of states that have not been sampled yet
            if (prevalence >= (1.0/numstates)){
              // update number of frames that contributed if the energy of the endstate was bellow its free energy difference + KbT
              DEBUG(1, "state 2 " << (eds_vi[state] - conf.current().energies.eds_vr) << "  Free ener "  << sim.param().eds.statefren[state]);
              if ((eds_vi[state] - conf.current().energies.eds_vr)  <= (sim.param().eds.statefren[state] + (1/beta))){
                  sim.param().eds.framecounts[state] += 1;
              }
              // check if we have seen enough frames for this state (convergence criteria) to update offsets
              //if (sim.param().eds.framecounts[state] > sim.param().eds.cc) {
              sim.param().eds.eir[state] -= prevalence/500;
              //}
            }
          } // loop over states
          for (unsigned int state = 1; state < numstates; state++) {
            sim.param().eds.eir[state] -= sim.param().eds.eir[0];
          }
          sim.param().eds.eir[0] = 0;
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

      // check if simulation should finish earlier due to convergence of adaptive search reached
      if (sim.param().eds.form == simulation::aeds_advanced_search || sim.param().eds.form == simulation::aeds_advanced_search2){
        bool convergence = true;
        for (unsigned int state = 0; state < numstates; state++) {
          if (sim.param().eds.framecounts[state] < sim.param().eds.cc){
            convergence = false;
            break;
          }
        }
        if (convergence){
          // if convergence is reached finish simulation
          return E_AEDS_CONVERGENCE;
        }

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


