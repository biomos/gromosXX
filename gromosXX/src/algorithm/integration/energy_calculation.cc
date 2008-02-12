/**
 * @file energy_calculation.cc
 * contains the implementation
 * for the Energy_Calculation class
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include "energy_calculation.h"

#include <util/error.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * calculate total energies
 * update the energy averages
 */
int algorithm::Energy_Calculation::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim,
 std::ostream & os,
 bool quiet
 )
{
  // os << "Energy calculation\n";
  
  // the resizing of the energy-arrays could be moved here...
  return 0;
}

/**
 * energy calculation
 */
int algorithm::Energy_Calculation::apply
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim
 )
{

  // update the energies
  if (conf.old().energies.calculate_totals()){
    std::cout << "\nError during MD run : energies are NaN!\n"
	      << std::endl;
    return E_NAN;
  }
  
  // perturbed energy derivatives
  if (sim.param().perturbation.perturbation){
    if (conf.old().perturbed_energy_derivatives.calculate_totals()){
      std::cout << "\nError during MD run : energy derivatives are NaN!\n"
		<< std::endl;
      return E_NAN;
    }

    conf.old().perturbed_energy_derivatives.entropy_term =
      conf.old().perturbed_energy_derivatives.potential_total *
      conf.old().energies.potential_total;
  }

  conf.current().averages.apply(topo, conf, sim);

  return 0;
}

