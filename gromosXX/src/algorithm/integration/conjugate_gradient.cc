/**
 * @file conjugate_gradient.cc
 * contains the implementation
 * for conjugate gradient energy minimisation
 * taken from GROMOS96 Manual
 * 
 * TODO:
 * 1. Test on different molecular systems
 * 2. Test posres with shake
 * 3. Revise changes in configuration.cc
 * 4. Unlock and test relevant feature pairs
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

#include "../../algorithm/algorithm/algorithm_sequence.h"
#include "../../interaction/forcefield/forcefield.h"
#include "../../algorithm/constraints/shake.h"
#include "../../algorithm/constraints/position_constraints.h"

#include "conjugate_gradient.h"

#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

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
    if (sim.param().minimise.ntem == 2)
      os << "\tFletcher-Reeves conjugate gradient\n";
    else if (sim.param().minimise.ntem == 3)
      os << "\tPolak-Ribiere conjugate gradient\n";

    os << "\tresetting search direction every n-th step : " << sim.param().minimise.ncyc;
    if (sim.param().minimise.ncyc == 0)
      os << " (no forced reset)";
    
    os << "\n" << std::scientific << std::setprecision(2)
       << "\trequested rms force in minimum             : " << sim.param().minimise.dele << "\n"
       << "\tinitial step size                          : " << sim.param().minimise.dx0 << "\n"
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
    if (sim.param().pairlist.skip_step > 1) {
    io::messages.add("For tight convergence, the pairlist should be generated every step or without cut-off",
            "Algorithm::conjugate_gradient",
            io::message::warning);
    }
  }
  // Set initial step size
  sim.minimisation_step_size() = sim.param().minimise.dx0;
  // Zero velocities
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;

  // Get pointers to algorithms
  // Position restraints
  do_posres = (sim.param().posrest.posrest == 3);
  if (do_posres) {
    cgrad_posres = dynamic_cast<algorithm::Position_Constraints *>(cgrad_seq.algorithm("Position_Constraints"));
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

  // SHAKE algorithm
  do_shake = (
    (sim.param().constraint.solute.algorithm == simulation::constr_shake)
    || (sim.param().system.nsm && sim.param().constraint.solvent.algorithm == simulation::constr_shake)
  );
  if (do_shake) {
    cgrad_shake = dynamic_cast<algorithm::Shake *>(cgrad_seq.algorithm("Shake"));
    if (cgrad_shake == NULL) { 
      std::cout << "Conjugate Gradient: could not get SHAKE algorithm"
                << "\n\t(internal error)" << std::endl;
      return 1;
    }
    // Initialise separate configuration for SHAKE
    conf_sh_init(conf);
  }

  // Initialise squared magnitude of search direction
  p_squared = 0.0;
  // Total number of iterations
  total_iterations = 0;

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
  #ifndef NDEBUG
    std::cout << std::scientific << std::setprecision(12);
  #else
    std::cout << std::scientific << std::setprecision(4);
  #endif
  DEBUG(10, "Step no. " << sim.steps());

  // Counter of interaction calculations
  unsigned counter_iter = 0;
  // Calculate the energies, as we need them to estimate a minimum along the search direction
  conf.current().energies.calculate_totals();
  
  // Step size for this step
  double step_size = sim.minimisation_step_size();
  DEBUG(10, "Current step size = " << step_size);

  // SHAKE forces
  if (do_shake) {
    double shake_step = Smin * sqrt(p_squared);
    if (shake_step == 0.0) {
      DEBUG(10,"Initial SHAKE of forces");
      shake_step = step_size;
    }
    if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
      io::messages.add("Unable to SHAKE forces of the lower boundary configuration", "Conjugate_Gradient", io::message::error);
      return error;
    }
  }

  // only if minimum number of steps were made
  if (sim.steps() > unsigned(sim.param().minimise.nmin)) {
    // Evaluate if the minimum RMS force criterion is met
    double f = 0.0;
    double f_max = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f, f_max);
    }
    f = sqrt(f / topo.num_atoms());
    f_max = sqrt(f_max);
    DEBUG(7, "RMS force = " << f << ", MAX force = " << f_max);
    DEBUG(7, "Total energy = " << conf.current().energies.potential_total + conf.current().energies.special_total);
    if (f < sim.param().minimise.dele) {
      std::cout << "CONJUGATE GRADIENT:\tMINIMUM REACHED\n";
      std::cout << "After " << total_iterations + 1 << " iterations \n";
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

  // SHAKE search directions
  if (do_shake && sim.steps()) {
    DEBUG(15, "Obtaining constrained search directions");
    shake_cgrads(topo, conf);
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

  // Calculate the search direction and the gradient
  double gradA = calculate_cgrad(topo, conf, beta);
  DEBUG(15, "gradA = " << gradA);

  // If gradA < 0.0, then A is beyond the minimum in the search direction,
  // so reset the search direction and calculate a new one
  if (gradA < 0.0) {
    DEBUG(1, "gradA below zero. Resetting the search direction");
    beta = 0.0;
    gradA = calculate_cgrad(topo, conf, beta);
    DEBUG(10, "After reset, gradA = " << gradA << ", beta = " << beta);
  }

  // Calculate upper boundary in the search direction
  // <p|p> Squared magnitude of search directions
  p_squared = 0.0;

  // If beta == 0, then p = f and <p|p> = gradA 
  if (beta == 0.0) {
    p_squared = gradA;
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      p_squared += math::abs2(conf.current().cgrad(i));
    }
  }
  DEBUG(15, "p_squared = " << p_squared);

  // Initial step size (chosen so rA->rB = step_size)
  double b_init = step_size / sqrt(p_squared);
  
  // Upper boundary for cubic interpolation
  double b = b_init;
  DEBUG(15, "b_init = " << b_init);

  // The energy of lower boundary equals the energy of initial configuration
  double eneA = conf.current().energies.potential_total + conf.current().energies.special_total;
  DEBUG(15, "eneA = " << eneA);

  // Current configuration becomes old
  conf.exchange_state();
  
  // Store intermediate interpolated configurations just for the interpolation criterion evaluation
  math::VArray posX_cur = conf.current().pos;
  math::VArray posX_old = conf.old().pos;
  
  // Gradient and energy in the upper boundary
  double gradB, eneB;
  
  // Counter of interval doubling
  unsigned counter_doub = 0;
  
  while(true) {
    // Calculate coordinates of the upper boundary configuration
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().pos(i) = conf.old().pos(i) + b * conf.old().cgrad(i);
    }

    // Perform interaction calculation
    if ( int error = evaluate_configuration(topo, conf, sim, eneB, false) ) {
      io::messages.add("Unable to evaluate the upper boundary configuration", "Conjugate_Gradient", io::message::error);
      return error;
    }

    counter_iter += 1;
    DEBUG(15, "eneB = " << eneB);

    if (do_shake) {
      double shake_step = step_size;
      if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
        io::messages.add("Unable to SHAKE forces of the upper boundary configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
    }

    // Calculate gradient along the search direction in the upper boundary as <p|f>
    gradB = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gradB += math::dot(conf.old().cgrad(i), conf.current().force(i));
    }
    DEBUG(15, "gradB = " << gradB);
    // If gradB < 0 or eneB > eneA, accept B as upper boundary for estimation of X
    if (gradB < 0 || eneB > eneA) {
      DEBUG(10, "Upper boundary accepted.");
      DEBUG(10, "b = " << b);
      break;
    }
    else {
      // Minimum is probably beyond B
      DEBUG(1, "Minimum is beyond upper boundary. Doubling the interval size.");
      b *= 2;
      DEBUG(10, "Increasing size of next step by 10%");
      sim.minimisation_step_size() *= 1.1;
      counter_doub += 1;
    }
  }
  double a = 0.0;
  double Z, W;
  double X, gradX, eneX;
  // Counter of interpolations
  unsigned counter_ipol = 0;
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
    // Calculate coordinates of the interpolated configuration
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().pos(i) = posX_cur(i) = conf.old().pos(i) + X * conf.old().cgrad(i);
    }

    // Calculate the X RMS displacement
    double disp_X = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      disp_X += math::abs2(posX_cur(i) - posX_old(i));
    }
    disp_X = sqrt(disp_X / topo.num_atoms());
    DEBUG(10, "RMS displacement of X = " << disp_X);

    // Accept X as new configuration, if criterion is met
    if (disp_X < sim.param().minimise.cgic || counter_ipol == unsigned(sim.param().minimise.cgim)) {
      //TEST
      /*if ( int error = evaluate_configuration(topo, conf, sim, eneX, do_shake) ) {
        io::messages.add("Unable to evaluate the interpolated configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
      // SHAKE forces of the interpolated configuration
      if (do_shake) {
        double shake_step = X * sqrt(p_squared);
        if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
          io::messages.add("Unable to SHAKE forces of the interpolated configuration", "Conjugate_Gradient", io::message::error);
          return error;
        }
      }

      // Calculate gradient along the search direction in X
      gradX = 0.0;
      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        gradX += math::dot(conf.old().cgrad(i), conf.current().force(i));
      }
      
      DEBUG(10, "Interpolation accepted\n"
        << "A = " << a << "\t eneA = " << eneA << "\t gradA = " << gradA << "\n"
        << "X = " << X << "\t eneX = " << eneX << "\t gradX = " << gradX << "\n"
        << "B = " << b << "\t eneB = " << eneB << "\t gradB = " << gradB
      );*/
      //END TEST
      break;
    }
    else {
      std::swap(posX_cur, posX_old);

      // Calculate interactions of the interpolated configuration
      if ( int error = evaluate_configuration(topo, conf, sim, eneX, false) ) {
        io::messages.add("Unable to evaluate the interpolated configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
      counter_iter += 1;
      DEBUG(10, "eneX = " << eneX);
      
      // SHAKE forces of the interpolated configuration
      if (do_shake) {
        double shake_step = X * sqrt(p_squared);
        if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
          io::messages.add("Unable to SHAKE forces of the interpolated configuration", "Conjugate_Gradient", io::message::error);
          return error;
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
      DEBUG(2, "New search interval: X-B");
          a = X;
          eneA = eneX;
          gradA = gradX;
      }
      // Otherwise move B to X
      else {
      DEBUG(2, "New search interval: A-X");
          b = X;
          eneB = eneX;
          gradB = gradX;
      }
    }
  }
  total_iterations += counter_iter + 1;
  Smin = X;
  DEBUG(7, "Minimum along the search direction accepted, Smin = " << Smin << "\n"
        << "Number of interval doubling: " << counter_doub << "\n"
        << "Number of cubic interpolations: " << counter_ipol << "\n"
        << "Number of interaction calculations: " << counter_iter << "\n"
        << "Total number of iterations so far: " << total_iterations);
        
  // If this is the last step, also perform one more calculation to print correct minimized energies
  if (sim.steps() == (unsigned(sim.param().step.number_of_steps) - 1)) {
    if ( int error = evaluate_configuration(topo, conf, sim, eneX, do_shake) ) {
      io::messages.add("Unable to evaluate the interpolated configuration", "Conjugate_Gradient", io::message::error);
      return error;
    }
    counter_iter += 1;
    if (do_shake) {
      double shake_step = X * sqrt(p_squared);
      if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
        io::messages.add("Unable to SHAKE forces of the interpolated configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
    }
    // Also print final RMS force
    double f = 0.0, f_max = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f, f_max);
    }
    f = sqrt(f / topo.num_atoms());
    f_max = sqrt(f_max);
    std::cout << "CONJUGATE GRADIENT:\tMIMIMUM CRITERION NOT MET\n";
    std::cout << "After " << total_iterations + 1 << " iterations \n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Final RMS force:\t" << f << "\n";
    std::cout << "Final MAX force:\t" << f_max << "\n";
  }
  else {
    // Modify the step size
    if (X < b_init / 10.0) {
      sim.minimisation_step_size() *= 0.9;

      DEBUG(15, "X below 10% of the search interval. Decreasing the step size by 10 percent.");
    }
    if (counter_doub && (sim.minimisation_step_size() > sim.param().minimise.dxm)) {
      sim.minimisation_step_size() = sim.param().minimise.dxm;
    }
  }
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  return 0;
}

/**
 * Shake current forces
 */
int algorithm::Conjugate_Gradient
::shake_forces
(
  topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  double shake_step
)
{
  DEBUG(7,"Obtaining constrained forces using SHAKE");
  // We use a separate configuration to not mess up the main configuration
  math::VArray & pos_old = conf_sh.old().pos;
  math::VArray & pos_cur = conf_sh.current().pos;
  // Forces of the main configuration
  math::VArray & force_cur = conf.current().force;

  // Copy current position
  pos_old = conf.current().pos;

  // <f|f> Squared magnitude of forces
  double f_squared = 0.0;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    f_squared += math::abs2(force_cur(i));
  }
  
  // Step size along unconstrained forces to SHAKE
  double step_sh = shake_step / sqrt(f_squared);

  DEBUG(10,"shake_forces: Total displacement = " << shake_step);
  DEBUG(15,"shake_forces: \n"
    << "f_squared = \t" << f_squared << "\n"
    << "step_sh = \t" << step_sh
  );
  
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    pos_cur(i) = pos_old(i) + step_sh * force_cur(i);
  }

  // Apply SHAKE
  if ( int error = cgrad_shake->apply(topo, conf_sh, sim) ) {
    io::messages.add("Unable to apply SHAKE algorithm", "Conjugate_Gradient", io::message::error);
    return error;
  }
  // Obtain and write constrained forces
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    force_cur(i) = (conf_sh.current().pos(i) - conf_sh.old().pos(i)) / step_sh;
  }
  return 0;
}

/**
 * Calculate the search direction coefficient
 */
double algorithm::Conjugate_Gradient
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
  else if (sim.param().minimise.ntem == 3) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f2 += math::dot(force_current(i) - force_old(i), force_current(i));
    }
  }
  DEBUG(15, "calculate_beta:" << "\n"
    << "f1 = " << f1 << "\n"
    << "f2 = " << f2 << "\n"
  );
  return f2 / f1;
}

/**
 * Update search direction and return the gradient
 */
double algorithm::Conjugate_Gradient
::calculate_cgrad
(
  const topology::Topology & topo,
  configuration::Configuration & conf,
  const double & beta
)
{
  const math::VArray & force = conf.current().force;
  const math::VArray & cgrad_old = conf.old().cgrad;
  math::VArray & cgrad_cur = conf.current().cgrad;
  double grad = 0.0;
  
  if (beta == 0.0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      cgrad_cur(i) = force(i);
      grad += math::abs2(force(i));
    }
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      cgrad_cur(i) = force(i) + beta * cgrad_old(i);
      grad += math::dot(cgrad_cur(i), force(i));
    }
  }
  return grad;
}

/**
 * Calculate interactions and energies of the conformation
 * Optionally also apply constraints
 */
int algorithm::Conjugate_Gradient
::evaluate_configuration
(
  topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  double & ene,
  bool do_pos_shake
)
{
  // If positional SHAKE is requested
  if (do_pos_shake) {
    DEBUG(7,"Obtaining constrained positions using SHAKE");
    if ( int error = cgrad_shake->apply(topo, conf, sim) ) {
      io::messages.add("Unable to apply SHAKE algorithm", "Conjugate_Gradient", io::message::error);
      return error;
    }
  }

  DEBUG(7,"Calculating interactions");
  if ( int error = cgrad_ff->apply(topo, conf, sim) ) {
      io::messages.add("Unable to calculate interactions", "Conjugate_Gradient", io::message::error);
    return error;
  }
  conf.current().energies.calculate_totals();
  ene = conf.current().energies.potential_total + conf.current().energies.special_total;

  if (do_posres) {
    DEBUG(7,"Applying position restraints");
    if ( int error = cgrad_posres->apply(topo, conf, sim) ) {
      io::messages.add("Unable to apply positional restraints", "Conjugate_Gradient", io::message::error);
      return error;
    }
  }

  return 0;
}

  
