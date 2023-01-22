/**
 * @file conjugate_gradient.cc
 * contains the implementation
 * for conjugate gradient energy minimisation
 * taken from GROMOS96 Manual
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
    io::messages.add("For tight convergence, the pairlist should be generated every step",
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
      std::cerr << "Conjugate Gradient: Could not get Position Constraints algorithm"
                << "\n\t(internal error)" << std::endl;
      return 1;
    }
  }

  // Forcefield to evaluate forces and energies
  cgrad_ff = dynamic_cast<interaction::Forcefield *>(cgrad_seq.algorithm("Forcefield"));
  if (cgrad_ff == NULL) {
    std::cerr << "Conjugate Gradient: could not get Interaction Calculation algorithm"
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
      std::cerr << "Conjugate Gradient: could not get SHAKE algorithm"
                << "\n\t(internal error)" << std::endl;
      return 1;
    }
    // Initialise separate configuration for SHAKE
    conf_sh_init(conf);
  }

  // Initialise squared magnitude of search direction
  p_squared = 0.0;
  // Counters of minimisation performance
  total_doubling = 0;
  total_interpolations = 0;
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
      return terminate(f, f_max, true);
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
  double beta = 0.0;
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

  // Calculate the search direction and the slope
  double gA = calculate_cgrad(topo, conf, beta);
  DEBUG(15, "gA = " << gA);

  // If gA < 0.0, then A is beyond the minimum in the search direction,
  // so reset the search direction and calculate a new one
  if (gA < 0.0) {
    DEBUG(1, "gA below zero. Resetting the search direction");
    beta = 0.0;
    gA = calculate_cgrad(topo, conf, beta);
    DEBUG(10, "After reset, gA = " << gA << ", beta = " << beta);
  }

  // Calculate upper boundary in the search direction
  // <p|p> Squared magnitude of search directions
  p_squared = 0.0;

  // If beta == 0, then p = f and <p|p> = gA 
  if (beta == 0.0) {
    if (gA == 0.0) {
      io::messages.add("Force is zero, unable to establish search direction", "Conjugate_Gradient", io::message::error);
      return E_NAN;
    }
    else if (math::gisnan(gA)) {
      io::messages.add("Force is NaN", "Conjugate_Gradient", io::message::error);
      return E_NAN;
    }
    p_squared = gA;
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      p_squared += math::abs2(conf.current().cgrad(i));
    }
  }
  DEBUG(15, "p_squared = " << p_squared);

  // Initial step size (chosen so |rB - rA| = step_size)
  double b_init = step_size / sqrt(p_squared);
  
  // Upper boundary for cubic interpolation
  double b = b_init;
  DEBUG(15, "b_init = " << b_init);

  // The energy of lower boundary equals the energy of initial configuration
  double eneA = conf.current().energies.potential_total + conf.current().energies.special_total;
  DEBUG(15, "eneA = " << eneA);

  // Current configuration becomes old
  conf.exchange_state();
  
  // Slope and energy in the upper boundary
  double gB = 0.0, eneB = 0.0;
  
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
      if ( int error = shake_forces(topo, conf, sim, step_size) ) {
        io::messages.add("Unable to SHAKE forces of the upper boundary configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
    }

    // Calculate slope along the search direction in the upper boundary as <p|f>
    gB = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gB += math::dot(conf.old().cgrad(i), conf.current().force(i));
    }

    DEBUG(15, "gB = " << gB);
    // If gB < 0 or eneB > eneA, accept B as upper boundary for estimation of X
    if (gB < 0 || eneB > eneA) {
      DEBUG(10, "Upper boundary accepted.");
      DEBUG(10, "b = " << b);
      break;
    }
    else {
      // Minimum is probably beyond B
      DEBUG(1, "Minimum is beyond upper boundary. Doubling the interval size.");
      b *= 2;
      DEBUG(10, "Increasing the size of next step by 10%");
      sim.minimisation_step_size() *= 1.1;
      counter_doub += 1;
    }
  }
  double a = 0.0;
  double Z = 0.0, W = 0.0;
  double X = 0.0;
  double gX = 0.0, eneX = 0.0, old_X = 0.0;
  // Counter of interpolations
  unsigned counter_ipol = 0;
  while(true) {
    old_X = X;
    // Estimation of X by cubic interpolation
    Z = (3 * (eneA - eneB) / (b - a)) - gA - gB;
    W = sqrt(Z * Z - gA * gB);
    X = b - (W - Z - gB) * (b - a) / (gA - gB + 2 * W);
    counter_ipol += 1;
    DEBUG(15, "X = " << X);
    if(math::gisnan(X)) { 
      io::messages.add("new coordinates are NaN", "Conjugate_Gradient", io::message::error);
    return E_NAN;
    }
    // Calculate coordinates of the interpolated configuration
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().pos(i) = conf.old().pos(i) + X * conf.old().cgrad(i);
    }

    // Calculate the X RMS displacement
    double disp_X = fabs(X - old_X) * sqrt(p_squared / topo.num_atoms());
    DEBUG(10, "RMS displacement of X = " << disp_X);

    // Accept X as new configuration, if criterion is met
    if (disp_X < sim.param().minimise.cgic || counter_ipol == unsigned(sim.param().minimise.cgim)) {
      break;
    }
    else {
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

      // Calculate slope along the search direction in X
      gX = 0.0;
      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        gX += math::dot(conf.old().cgrad(i), conf.current().force(i));
      }
      DEBUG(10, "gX = " << gX);
      
      // If energy descends in X and eneX < eneA, move A to X
      if (gX > 0 && eneX < eneA) {
      DEBUG(2, "New search interval: X-B");
          a = X;
          eneA = eneX;
          gA = gX;
      }
      // Otherwise move B to X
      else {
      DEBUG(2, "New search interval: A-X");
          b = X;
          eneB = eneX;
          gB = gX;
      }
    }
  }
  Smin = X;

  total_doubling += counter_doub;
  total_interpolations += counter_ipol;
  total_iterations += counter_iter + 1;

  DEBUG(7, "Minimum along the search direction accepted, Smin = " << Smin);
  DEBUG(15,
           "Number of interval doubling: " << counter_doub << "\n"
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
    // TEST
    DEBUG(7, "Total energy = " << conf.current().energies.potential_total + conf.current().energies.special_total)
    DEBUG(7, "El (RF) = " << conf.current().energies.crf_total)
    // END TEST

    if (do_shake) {
      double shake_step = X * sqrt(p_squared);
      if ( int error = shake_forces(topo, conf, sim, shake_step) ) {
        io::messages.add("Unable to SHAKE forces of the interpolated configuration", "Conjugate_Gradient", io::message::error);
        return error;
      }
    }
    // Also print the final RMS force
    double f = 0.0, f_max = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f, f_max);
    }
    f = sqrt(f / topo.num_atoms());
    f_max = sqrt(f_max);
    // Terminate
    return terminate(f, f_max, false);
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
  
  // Normalize step size along unconstrained forces
  double shake_step_norm = shake_step / sqrt(f_squared);
  DEBUG(10,"shake_forces: Total displacement = " << shake_step);
  DEBUG(15,"shake_forces: \n"
    << "f_squared = \t" << f_squared << "\n"
    << "shake_step_norm = \t" << shake_step_norm
  );
  
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    pos_cur(i) = pos_old(i) + shake_step_norm * force_cur(i);
  }

  // Apply SHAKE
  if ( int error = cgrad_shake->apply(topo, conf_sh, sim) ) {
    io::messages.add("Unable to apply SHAKE algorithm", "Conjugate_Gradient", io::message::error);
    return error;
  }
  // Obtain and write constrained forces
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    force_cur(i) = (conf_sh.current().pos(i) - conf_sh.old().pos(i)) / shake_step_norm;
  }
  return 0;
}

/**
 * Calculate the search direction coefficient beta
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
 * Update search direction and return the slope in the search dierction
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
  double slope = 0.0;
  
  if (beta == 0.0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      cgrad_cur(i) = force(i);
      slope += math::abs2(force(i));
    }
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      cgrad_cur(i) = force(i) + beta * cgrad_old(i);
      slope += math::dot(cgrad_cur(i), force(i));
    }
  }
  return slope;
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

/**
 * Terminate the minimisation
 */
/** This could be done later with MINIMISATION block */
int algorithm::Conjugate_Gradient
::terminate
(
  double rms_force,
  double max_force,
  int minimum
)
{
  int error = 0;
  std::cout << std::string(60, '-') << "\n";
  if (minimum) {
    error = E_MINIMUM_REACHED;
    std::cout << "CONJUGATE GRADIENT:\tMINIMUM REACHED\n";
  }
  else {
    std::cout << "CONJUGATE GRADIENT:\tMIMIMUM CRITERION NOT MET\n";
  }
  std::cout << "Total interaction calculations : " << std::setw(10) << std::right << total_iterations + 1 << "\n";
  std::cout << "Total search interval doublings: " << std::setw(10) << std::right << total_doubling << "\n";
  std::cout << "Total interpolations           : " << std::setw(10) << std::right << total_interpolations << "\n";
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "Final RMS force                : " << rms_force << "\n";
  std::cout << "Final MAX force                : " << max_force << "\n";
  std::cout << std::string(60, '-') << "\n\n";
  return error;
}

  
