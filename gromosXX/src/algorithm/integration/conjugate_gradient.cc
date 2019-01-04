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
 * 4. Test on large number of different molecular systems
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

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
    if (sim.param().minimise.ntem == 2) {
      os << "\tFletcher-Reeves conjugate gradient\n";
    }
    else {
      os << "\tPolak-Ribiere conjugate gradient\n";
    }
    os << "\tresetting search direction every n-th step : " << sim.param().minimise.ncyc << "\n"
       << std::scientific
       << "\trequested rms force in minimum             : " << sim.param().minimise.dele << "\n" 
       << std::fixed
       << "\tminimum and starting step size             : " << sim.param().minimise.dx0 << "\n"
       << "\tmaximum step size                          : " << sim.param().minimise.dxm << "\n"
       << "\tminimum steps                              : " << sim.param().minimise.nmin << "\n"
       << "\tmaximum cubic interpolations per step      : " << sim.param().minimise.cgim << "\n"
       << std::scientific
       << "\tdisplacement criterion on interpolation    : " << sim.param().minimise.cgic << "\n"
       << std::fixed
       ;
  }
  if (!quiet) {
    if (sim.param().minimise.flim != 0) {
      os  << "\tlimiting the force to    : "
	        << sim.param().minimise.flim << "\n";
    }
    os << "END\n";
    if (sim.param().pairlist.skip_step > 1 && sim.param().pairlist.skip_step < sim.param().step.number_of_steps) {
    io::messages.add("For tight convergence, the pairlist should be generated every step",
            "Algorithm::conjugate_gradient",
            io::message::warning);
    }
  }
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
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
  int unsigned its = 0; // Counter of interaction calculations
  DEBUG(15,"step no.:\t" << sim.steps());
  // Calculate the energies, as we need them to estimate a minimum along the search direction
  conf.current().energies.calculate_totals();
  // Keep the step size above user defined size
  if (sim.minimisation_step_size() < sim.param().minimise.dx0) {
    sim.minimisation_step_size() = sim.param().minimise.dx0;
  }
  DEBUG(10, "Current step size = " << sim.minimisation_step_size()
);
  // only if minimum number of steps were made
  if (sim.steps() > unsigned(sim.param().minimise.nmin)) {
    // check whether minimum is reached by the RMS force criterion
    double f = 0.0, f_max = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f,f_max);
    }
    f = sqrt(f/topo.num_atoms());
    f_max = sqrt(f_max);
    DEBUG(10, "RMS force = " << f << ", MAX force = " << f_max);
    DEBUG(15, "Total energy = " << conf.current().energies.potential_total + conf.current().energies.special_total);
    if (f < sim.param().minimise.dele) {
      std::cout << "CONJUGATE GRADIENT:\tMINIMUM REACHED\n";
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
// Check, whether we have any non-zero old force, otherwise we dont have to calculate beta
  bool no_force = true;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    if (math::abs2(conf.old().force(i)) > math::epsilon) {
      no_force = false;
      break;
    }
  }
  // If old forces are non-zero and this is not a resetting step, calculate beta
  double beta = 0.0;
  if (!no_force && (sim.param().minimise.ncyc == 0 || sim.steps() % sim.param().minimise.ncyc != 0)) {
    beta = calculate_beta(topo,conf,sim);
    DEBUG(15, "beta = " << beta);
  }
  // Otherwise reset the search direction by keeping beta = 0
  else {
      DEBUG(1,
      "(Re)initializing the conjugate gradient search direction\n"
      << "beta = " << beta);
  }
  // Calculate the search direction
  double b = calculate_cgrad(topo,conf,beta);
  // Calculate a gradient of energy along the search direction
  double gA;
  // If beta = 0, gA is identical to b
  if (beta == 0.0) {
    gA = b;
  }
  else {
    gA = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gA += math::dot(conf.current().cgrad(i), conf.current().force(i));
    }
  }
  DEBUG(15, "gA = " << gA);
  // If gA < 0.0, then A is beyond the minimum in the search direction,
  // so reset the search direction and calculate a new one
  if (gA < 0.0) {
    DEBUG(1, "gA below zero. Resetting the search direction");
    beta = 0.0;
    b = calculate_cgrad(topo,conf,beta);
    gA = b;
    DEBUG(10, "After reset, gA = " << gA << ", beta = " << beta);
  }
  // Calculate upper boundary in the search direction
  double init_b = b = sim.minimisation_step_size() / sqrt(b);
  DEBUG(10, "b = " << b);
  // Energy of lower boundary equals energy of initial conf
  double eneA = conf.current().energies.potential_total + conf.current().energies.special_total;
  DEBUG(15, "eneA = " << eneA);
  // Create confX to store intermediate configurations and evaluate interpolations
  configuration::Configuration confX = conf;
  double gB, eneB, X;
  unsigned int doubled = 0; // Counter of interval doubling
  unsigned int ints = 0; // Number of interpolations done
  conf.exchange_state();
  while (true) {
    // Calculate new coordinates, energies, forces in upper boundary
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().pos(i) = conf.old().pos(i) + b * conf.old().cgrad(i);
    }
    cg_ff.apply(topo, conf, sim);    
    its += 1;
    conf.current().energies.calculate_totals();
    eneB = conf.current().energies.potential_total + conf.current().energies.special_total;
    DEBUG(15, "eneB = " << eneB);
    // Calculate gradient along the search direction in upper boundary
    gB = 0.0;
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      gB += math::dot(conf.old().cgrad(i), conf.current().force(i));
    }
    DEBUG(15, "gB = " << gB);
    // If gB < 0 or eneB > eneA, accept b as upper boundary for estimation of X
    if (gB < 0 || eneB > eneA) {
      double a = 0.0;
      double Z, W;
      double gX, eneX;
      while(true) {
        // Estimation of X by cubic interpolation
        Z = (3 * (eneA - eneB) / (b - a)) - gA - gB;
        W = sqrt(Z * Z - gA * gB);
        X = b - (W - Z - gB) * (b - a) / (gA - gB + 2 * W);
        ints += 1;
        DEBUG(15, "X = " << X);
        if(math::gisnan(X)) {
          io::messages.add("new coordinates are NaN", "Conjugate_Gradient", io::message::error);
        return E_NAN;
        }
        // Calculate new coordinates of X
        confX.exchange_state();
        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          conf.current().pos(i) = confX.current().pos(i) = conf.old().pos(i) + X * conf.old().cgrad(i);
        }     
        // Calculate RMS displacement of X
        double d = 0.0;
        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          d += math::abs2(confX.current().pos(i) - confX.old().pos(i));
        }
        d = sqrt(d/topo.num_atoms());
        DEBUG(10, "RMS displacement of X = " << d);
        // Accept X as new configuration, if criterion is met
        if (d < sim.param().minimise.cgic || ints == sim.param().minimise.cgim) {
          break;
        }
        // Calculate forces and energies of X
        cg_ff.apply(topo, conf, sim);
        its += 1;
        conf.current().energies.calculate_totals();
        eneX = conf.current().energies.potential_total + conf.current().energies.special_total;
        DEBUG(10, "eneX = " << eneX);
        // Calculate gradient along the search direction in X
        gX = 0.0;
        for(unsigned int i=0; i<topo.num_atoms(); ++i) {
          gX += math::dot(conf.old().cgrad(i), conf.current().force(i));
        }
        DEBUG(10, "gX = " << gX);
        // If energy descends in X and eneX < eneA, move A to X
        if (gX > 0 && eneX < eneA) {
        DEBUG(2, "Moving boundary A to X");
            a = X;
            eneA = eneX;
            gA = gX;
        // Otherwise move B to X
        } else {
        DEBUG(2, "Moving boundary B to X");
            b = X;
            eneB = eneX;
            gB = gX;
        }
      }
      break;
    }
    else { // Minimum is probably beyond b
      DEBUG(1, "Minimum is beyond upper boundary, Doubling the interval size");
      b *= 2;
      DEBUG(15, "Increasing the next step size by 10%");
      sim.minimisation_step_size() *= 1.1;
      doubled += 1;
    }
  }
  DEBUG(7, "Minimum along the search direction accepted, X = " << X << "\n"
  << "Times the search interval size was doubled: " << doubled << "\n"
  << "Number of cubic interpolations: " << ints << "\n"
  << "Total number of interactions calculation: " << its);
  // If this is the last step, also perform one more calculation to print correct minimized energies
  if (sim.steps() == (unsigned(sim.param().step.number_of_steps) - 1)) {
    cg_ff.apply(topo, conf, sim);
    its += 1;
    conf.current().energies.calculate_totals();
    double f = 0.0, f_max = 0.0;
    // Also print final RMS force
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f += math::abs2(conf.current().force(i));
      f_max = std::max(f,f_max);
    }
    f = sqrt(f/topo.num_atoms());
    f_max = sqrt(f_max);
    std::cout << "CONJUGATE GRADIENT:\tMIMIMUM CRITERION NOT MET\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Final RMS force:\t" << f << "\n";
    std::cout << "Final MAX force:\t" << f_max << "\n";
  } else {
    // Modify the step size
    if (X < init_b / 5.0) {
      sim.minimisation_step_size() *= 0.9;
      DEBUG(15, "X below one fifth of interval. Decreasing the step size by 10 percent.");
    }
    if (doubled && sim.minimisation_step_size() > sim.param().minimise.dxm) {
      sim.minimisation_step_size() = sim.param().minimise.dxm;
    }
  }
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  return 0;
}

inline double algorithm::Conjugate_Gradient
::calculate_beta
(
  topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim
)
{ 
  double f1 = 0.0, f2 = 0.0;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    f1 += math::abs2(conf.old().force(i));
  }
  // Fletcher-Reeves
  if (sim.param().minimise.ntem == 2) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f2 += math::abs2(conf.current().force(i));
    }
  }
  // Polak-Ribiere
  else { //sim.param().minimise.ntem == 3
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      f2 += math::dot(conf.current().force(i) - conf.old().force(i), conf.current().force(i));
    }
  }
  return f2 / f1;
}

inline double algorithm::Conjugate_Gradient
::calculate_cgrad
(
  topology::Topology & topo,
  configuration::Configuration & conf,
  double & beta
)
{
  double b = 0;
  if (beta == 0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        conf.current().cgrad(i) = conf.current().force(i);
        b += math::abs2(conf.current().force(i));
    }
  }
  else {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      conf.current().cgrad(i) = conf.current().force(i) + beta * conf.old().cgrad(i);
      b += math::abs2(conf.current().cgrad(i));
    }
  }
  return b;
}
