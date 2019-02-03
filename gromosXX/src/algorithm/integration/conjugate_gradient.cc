/**
 * @file conjugate_gradient.cc
 * contains the implementation
 * for conjugate gradient energy minimisation
 * taken from GROMOS96 Manual
 * 
 * TODO:
 * 1. Decide what to print out (as in GROMOS96?)
 * 2. Unlock supported and meaningful features
 * 3. Allow SHAKE
 * 4. Test on different molecular systems
 * 5. Pairlist generation issue - is warning enough?
 * 6. No need to parallelization - this code scales with ~N, while energy calculation scales with ~N^2
 * 7. Test posres with shake
 * 8. For efficiency include only SHAKEN solvent/solute in constrained force calculation
 * 
 * WE SHOULD DO POSRES AND SHAKE AFTER EVERY CHANGE OF POSITIONS
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

// TEST
#include "../../algorithm/algorithm/algorithm_sequence.h"
//#include "../io/topology/in_topology.h"
//#include "../../algorithm/constraints/create_constraints.h"
//#include "../../algorithm/constraints/create_constraints.h"
#include "../../math/periodicity.h"
#include "../../algorithm/constraints/shake.h"
#include "../../algorithm/constraints/position_constraints.h"
// END TEST
#include "conjugate_gradient.h"

#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

// Additional algorithms used in the step
/**
 * Positional constraints 
 */
algorithm::Position_Constraints * cgrad_posres;

/**
 * SHAKE algorithm
 */
algorithm::Shake * cgrad_shake;

/**
 * Calculation of interactions
 */
interaction::Forcefield * cgrad_ff;

/**
 * conjugate gradient step.
 */
int algorithm::Conjugate_Gradient
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet) {
    os << "ENERGY MINIMISATION\n";
    if (sim.param().minimise.ntem == 2) {
      os << "\tFletcher-Reeves conjugate gradient\n";
    }
    else os << "\tPolak-Ribiere conjugate gradient\n";

    os << "\tresetting search direction every n-th step : " << sim.param().minimise.ncyc;
    if (sim.param().minimise.ncyc == 0) {
      os << " (no forced reset)";
    }
    os << "\n" << std::scientific << std::setprecision(2)
       << "\trequested rms force in minimum             : " << sim.param().minimise.dele << "\n"
       << "\tminimum and starting step size             : " << sim.param().minimise.dx0 << "\n"
       << "\tmaximum step size                          : " << sim.param().minimise.dxm << "\n"
       << "\tminimum steps                              : " << sim.param().minimise.nmin << "\n"
       << "\tmaximum cubic interpolations per step      : " << sim.param().minimise.cgim << "\n"
       << "\tdisplacement criterion on interpolation    : " << sim.param().minimise.cgic << "\n"
       << std::fixed << std::setprecision(4);
       ;
  }
  if (!quiet) {
    if (sim.param().minimise.flim != 0) {
      os  << "\tlimiting the force to    : "
	        << sim.param().minimise.flim << "\n";
    }
    os << "END\n";
    if (sim.param().pairlist.skip_step > 1 && sim.param().pairlist.skip_step < sim.param().step.number_of_steps) {
    io::messages.add("For tight convergence, the pairlist should be generated every step or without cut-off",
            "Algorithm::conjugate_gradient",
            io::message::warning);
    }
  }
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;

  // Get pointers to algorithms and perform some basic checks
  // Position restraints
  do_posres = sim.param().posrest.posrest == 3;
  cgrad_posres = dynamic_cast<algorithm::Position_Constraints *>(cgrad_seq.algorithm("Position_Constraints"));
  if (do_posres) {
    if (cgrad_posres == NULL) {
      std::cout << "Conjugate Gradient: Could not get Position Constraints algorithm"
                << "\n\t(internal error)" << std::endl;
      return 1;
    }
  }
  // Forcefield to evaluate forces and energies
  cgrad_ff = dynamic_cast<interaction::Forcefield *>(cgrad_seq.algorithm("Forcefield"));
  if (cgrad_ff == NULL) {
    std::cout << "Conjugate Gradient: could not get Interaction Calculation algorithm"
              << "\n\t(internal error)" << std::endl;
    return 1;
  }
  // Shake algorithm
  do_shake = (
    sim.param().constraint.solute.algorithm == simulation::constr_shake
    || (sim.param().system.nsm && sim.param().constraint.solvent.algorithm == simulation::constr_shake)
  );
  cgrad_shake = dynamic_cast<algorithm::Shake *>(cgrad_seq.algorithm("Shake"));
  if (do_shake) {
    if (cgrad_shake == NULL) { 
      std::cout << "Conjugate Gradient: could not get SHAKE algorithm"
                << "\n\t(internal error)" << std::endl;
      return 1;
    }
  }
  return 0;
}

/**
 * conjugate gradient step.
 */
int algorithm::Conjugate_Gradient
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  //int error;
  #ifndef NDEBUG
    std::cout << std::scientific << std::setprecision(12);
  #else
    std::cout << std::scientific << std::setprecision(4);
  #endif
  // Counter of interaction calculations
  int counter_iter = 0;
  DEBUG(15,"step no.:\t" << sim.steps());
  // Calculate the energies, as we need them to estimate a minimum along the search direction
  conf.current().energies.calculate_totals();
  // Keep the step size above user defined size
  if (sim.minimisation_step_size() < sim.param().minimise.dx0) {
    sim.minimisation_step_size() = sim.param().minimise.dx0;
  }
  DEBUG(10, "Current step size = " << sim.minimisation_step_size());


  if (do_shake) {
    // Obtain initial constraint forces
    configuration::Configuration conf_sh = conf;
    conf_sh.exchange_state();
    double f_squared = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f_squared += math::abs2(conf_sh.old().force(i));
    }
    // Special initial step with SHAKE
    double step_sh = sim.minimisation_step_size() / sqrt(f_squared);
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf_sh.current().pos(i) = conf_sh.old().pos(i) + step_sh * conf_sh.old().force(i);
    }
    // TODO: ALSO APPLY POSRES???
    cgrad_shake->apply(topo, conf_sh, sim);
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().force(i) = (conf_sh.current().pos(i) - conf_sh.old().pos(i)) / step_sh;
    }
  }

  // only if minimum number of steps were made
  if (sim.steps() > unsigned(sim.param().minimise.nmin)) {
    // check whether minimum is reached by the RMS force criterion
    double f = 0.0, f_max = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f, f_max);
    }
    f = sqrt(f / topo.num_atoms());
    f_max = sqrt(f_max);
    DEBUG(10, "RMS force = " << f << ", MAX force = " << f_max);
    DEBUG(15, "Total energy = " << conf.current().energies.potential_total + conf.current().energies.special_total);
    if (f < sim.param().minimise.dele) {
      std::cout << "CONJUGATE GRADIENT:\tMINIMUM REACHED\n";
      std::cout << "After " << sim.minimisation_iterations() + 1 << " iterations \n";
      std::cout << "Final RMS force:\t" << f << "\n";
      std::cout << "Final MAX force:\t" << f_max << "\n";
      return E_MINIMUM_REACHED;
    }
  }
  // limit the maximum force
  if (sim.param().minimise.flim != 0.0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      const double fs = math::abs(conf.current().force(i));
      if (fs > sim.param().minimise.flim) {
        DEBUG(15,"fs = " << fs 
             << ", flim = " << sim.param().minimise.flim);
        DEBUG(15,"force (" << i << ") = " << math::v2s(conf.current().force(i))
             << ", factor = " << sim.param().minimise.flim / fs);
        conf.current().force(i) *= sim.param().minimise.flim / fs;
      }
    }
  }
// Check, whether we have any non-zero old force, otherwise we cannot calculate beta
  bool no_force = true;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    if (math::abs2(conf.old().force(i)) > 0.0) {
      no_force = false;
      break;
    }
  }
  // If old forces are non-zero and this is not a resetting step, calculate beta
  double beta;
  if (!no_force && (sim.param().minimise.ncyc == 0 || sim.steps() % sim.param().minimise.ncyc != 0)) {
    beta = calculate_beta(topo, conf, sim);
    DEBUG(15, "beta = " << beta);
  }
  // Otherwise reset the search direction by keeping beta = 0
  else {
    beta = 0.0;
    DEBUG(1,
      "(Re)initializing the conjugate gradient search direction\n"
      << "beta = " << beta);
  }
  // Calculate the search direction and <p|p>
  double p_squared = calculate_cgrad(topo, conf, beta);
  // Calculate a gradient along the search direction
  double gradA;
  // If beta = 0, gradA is identical to p_squared
  if (beta == 0.0) {
    gradA = p_squared;
  }
  else {
    gradA = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gradA += math::dot(conf.current().cgrad(i), conf.current().force(i));
    }
  }
  DEBUG(15, "gradA = " << gradA);
  // If gradA < 0.0, then A is beyond the minimum in the search direction,
  // so reset the search direction and calculate a new one
  if (gradA < 0.0) {
    DEBUG(1, "gradA below zero. Resetting the search direction");
    beta = 0.0;
    p_squared = calculate_cgrad(topo, conf, beta);
    gradA = p_squared;
    DEBUG(10, "After reset, gradA = " << gradA << ", beta = " << beta);
  }
  // Calculate upper boundary in the search direction
  double b_init = sim.minimisation_step_size() / sqrt(p_squared);
  double b = b_init;
  DEBUG(10, "b = " << b);
  // The energy of lower boundary equals the energy of initial conf
  double eneA = conf.current().energies.potential_total + conf.current().energies.special_total;
  DEBUG(15, "eneA = " << eneA);
  // Create confX to store intermediate configurations and evaluate interpolations
  //configuration::Configuration confB, confX = conf;
  configuration::Configuration confX = conf;
  double gradB, eneB;
  // Counter of interval doubling
  int counter_doub = 0;
  
  conf.exchange_state();

  //double deltaX_uc, abs_p, abs_f_uc;
  while(true) {

    // Calculate new coordinates in upper boundary
    // TODO: MAKE FUNCTION OF THIS update_pos(new_positions, const old_positions, const coefficient, const search_vector)
    
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().pos(i) = conf.old().pos(i) + b * conf.old().cgrad(i);
    }

    if (do_posres) {
        // posres just simply restores old positions of restrained atoms and zeroes forces
        DEBUG(15, "Calling Positional Restraints");
        cgrad_posres->apply(topo, conf, sim);
    }

    cgrad_ff->apply(topo,conf,sim);
    counter_iter += 1;

    if (do_shake) {
      configuration::Configuration conf_sh = conf;
      conf_sh.exchange_state();
      double step_sh, f_squared = 0.0;
      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        f_squared += math::abs2(conf_sh.old().force(i));
      }
      step_sh = b * sqrt(p_squared / f_squared);

      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        conf_sh.current().pos(i) = conf_sh.old().pos(i) + step_sh * conf_sh.old().force(i);
      }
      cgrad_shake->apply(topo,conf_sh,sim);
      //conf.current().pos = conf_sh.current().pos; // TODO: SHOULD WE UPDATE POS AT B?
        
      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        conf.current().force(i) = (conf_sh.current().pos(i) - conf_sh.old().pos(i)) / step_sh;
      }
    }
    
    conf.current().energies.calculate_totals();
    eneB = conf.current().energies.potential_total + conf.current().energies.special_total;
    DEBUG(15, "eneB = " << eneB);
    // Calculate gradient along the search direction in upper boundary as <p|f>
    gradB = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gradB += math::dot(conf.old().cgrad(i), conf.current().force(i));
    }
    DEBUG(15, "gradB = " << gradB);
    // If gradB < 0 or eneB > eneA, accept B as upper boundary for estimation of X
    if (gradB < 0 || eneB > eneA) { // TODO: ADD MAX ITERATIONS?
      DEBUG(15, "Upper boundary accepted.");
      break;
    }
    else { // Minimum is probably beyond B
      DEBUG(1, "Minimum is beyond upper boundary. Doubling the interval size.");
      b *= 2;
      DEBUG(15, "Increasing the next step size by 10%");
      sim.minimisation_step_size() *= 1.1;
      counter_doub += 1;
    }
  }
  double a = 0.0;
  double Z, W;
  double X, gradX, eneX;
  // Counter of interpolations
  int counter_ipol = 0;
  while(true) {
    // Estimation of X by cubic interpolation
    Z = (3 * (eneA - eneB) / (b - a)) - gradA - gradB;
    W = sqrt(Z * Z - gradA * gradB);
    X = b - (W - Z - gradB) * (b - a) / (gradA - gradB + 2 * W);
    counter_ipol += 1;
    DEBUG(15, "X = " << X);
    if(math::gisnan(X)) {
      io::messages.add("new coordinates are NaN", "Conjugate_Gradient", io::message::error);
    return E_NAN;
    }
    // Calculate new coordinates of X
    confX.exchange_state();
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      // confX is used only as reference for threshold
      conf.current().pos(i) = confX.current().pos(i) = conf.old().pos(i) + X * conf.old().cgrad(i);
    }
    

    if (do_posres) {
      // posres just simply restores old positions of restrained atoms and zeroes forces
      DEBUG(15, "Calling Positional Restraints");
      cgrad_posres->apply(topo, conf, sim);
    }
    if (do_shake) {
      cgrad_shake->apply(topo,conf,sim);
      confX.current().pos = conf.current().pos;
    }

    // Calculate the X RMS displacement
    double disp_X = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      disp_X += math::abs2(confX.current().pos(i) - confX.old().pos(i));
    }
    disp_X = sqrt(disp_X / topo.num_atoms());
    DEBUG(10, "RMS displacement of X = " << disp_X);
    // Accept X as new configuration, if criterion is met
    if (disp_X < sim.param().minimise.cgic || counter_ipol == sim.param().minimise.cgim) {
      break;
    }

    cgrad_ff->apply(topo,conf,sim);
    counter_iter += 1;
    conf.current().energies.calculate_totals();
    eneX = conf.current().energies.potential_total + conf.current().energies.special_total;
    DEBUG(10, "eneX = " << eneX);
    


    // Calculate the energy of new constrained positions
    {
      if (do_posres) {
          // posres just simply restores old positions of restrained atoms and zeroes forces
          DEBUG(15, "Calling Positional Restraints");
          cgrad_posres->apply(topo, conf, sim);
      }

      if (do_shake) {
        configuration::Configuration conf_sh = conf;
        conf_sh.exchange_state();
        double step_sh, f_squared = 0.0;
        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          f_squared += math::abs2(conf_sh.old().force(i));
        }
        step_sh = X * sqrt(p_squared / f_squared);

        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          conf_sh.current().pos(i) = conf_sh.old().pos(i) + step_sh * conf_sh.old().force(i);
        }
        cgrad_shake->apply(topo, conf_sh, sim);
        //conf.current().pos = conf_sh.current().pos; // TODO: SHOULD WE UPDATE POS AT B?
        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          conf.current().force(i) = (conf_sh.current().pos(i) - conf_sh.old().pos(i)) / step_sh;
        }
      }
    }
    // Calculate gradient along the search direction in X
    gradX = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gradX += math::dot(conf.old().cgrad(i), conf.current().force(i));
    }
    DEBUG(10, "gradX = " << gradX);
    // If energy descends in X and eneX < eneA, move A to X
    if (gradX > 0 && eneX < eneA) {
    DEBUG(2, "Moving boundary A to X");
        a = X;
        eneA = eneX;
        gradA = gradX;
    // Otherwise move B to X
    }
    else {
    DEBUG(2, "Moving boundary B to X");
        b = X;
        eneB = eneX;
        gradB = gradX;
    }
  }
  if (do_shake) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.old().cgrad(i) = (conf.current().pos(i) - conf.old().pos(i)) / X;
    }
  }
  sim.minimisation_iterations() += counter_iter + 1;
  DEBUG(7, "Minimum along the search direction accepted, X = " << X << "\n"
        << "Number of interval doubling: " << counter_doub << "\n"
        << "Number of cubic interpolations: " << counter_ipol << "\n"
        << "Total number of interaction calculations: " << counter_iter);
  // If this is the last step, also perform one more calculation to print correct minimized energies
  if (sim.steps() == (unsigned(sim.param().step.number_of_steps) - 1)) {
    cgrad_ff->apply(topo, conf, sim);
    //if (0 != (error = evaluate_conf(topo, conf, sim))) {
    //  return error;
    //}
    counter_iter += 1;
    conf.current().energies.calculate_totals();
    double f = 0.0, f_max = 0.0;
    // Also print final RMS force
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f, f_max);
    }
    f = sqrt(f / topo.num_atoms());
    f_max = sqrt(f_max);
    std::cout << "CONJUGATE GRADIENT:\tMIMIMUM CRITERION NOT MET\n";
    std::cout << "After " << sim.minimisation_iterations() + 1 << " iterations \n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Final RMS force:\t" << f << "\n";
    std::cout << "Final MAX force:\t" << f_max << "\n";
  } else {
    // Modify the step size
    if (X < b_init / 5.0) {
      sim.minimisation_step_size() *= 0.9;
      DEBUG(15, "X below one fifth of interval. Decreasing the step size by 10 percent.");
    }
    if (counter_doub && sim.minimisation_step_size() > sim.param().minimise.dxm) {
      sim.minimisation_step_size() = sim.param().minimise.dxm;
    }
  }
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  return 0;
}
/**
 * Calculate the search direction coefficient
 */
inline double algorithm::Conjugate_Gradient
::calculate_beta
(
  const topology::Topology & topo,
  const configuration::Configuration & conf,
  const simulation::Simulation & sim
)
{
  const math::VArray & force_old = conf.old().force;
  const math::VArray & force_current = conf.current().force;

  double f1 = 0.0, f2 = 0.0;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    f1 += math::abs2(force_old(i));
  }
  // Fletcher-Reeves
  if (sim.param().minimise.ntem == 2) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f2 += math::abs2(force_current(i));
    }
  }
  // Polak-Ribiere
  else { //sim.param().minimise.ntem == 3
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f2 += math::dot(force_current(i) - force_old(i), force_current(i));
    }
  }
  DEBUG(15, "f2 / f1 =" << f2 / f1);
  return f2 / f1;
}

/**
 * Update search directions and also return sum of their
 * squared sizes to be used for the search interval size
 */
inline double algorithm::Conjugate_Gradient
::calculate_cgrad
(
  const topology::Topology & topo,
  configuration::Configuration & conf,
  const double & beta
)
{
  const math::VArray & force = conf.current().force;
  double p_squared = 0.0;
  if (beta == 0.0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        conf.current().cgrad(i) = force(i);
        p_squared += math::abs2(force(i));
    }
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().cgrad(i) = force(i) + beta * conf.old().cgrad(i);
      p_squared += math::abs2(conf.current().cgrad(i));
    }
  }
  return p_squared;
}
/**
 * Calculate interactions and energies of the conformation
 * Optionally also apply constraints
 */
inline int algorithm::Conjugate_Gradient
::evaluate_conf
(
  topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim
)
{
  DEBUG(15,"Evaluating configuration");
  int error;
  if (do_posres && 0 != (error = cgrad_posres->apply(topo, conf, sim))) {
    std::cout << "Conjugate Gradient: Error applying positional constraints" << std::endl;
    return error;
  };
  if (do_shake && 0 != (error = cgrad_shake->apply(topo, conf, sim))) {
    std::cout << "Conjugate Gradient: Error applying SHAKE algorithm" << std::endl;
    return error;
  };
  if (0 != (error = cgrad_ff->apply(topo, conf, sim))) {
    std::cout << "Conjugate Gradient: Error calculating interactions" << std::endl;
    return error;
  };
  return 0;
}

  
