/**
 * @file check_forcefield.cc
 */


#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <simulation/parameter.h>
#include <interaction/interaction_types.h>
#include <io/instream.h>
#include <util/parse_tcouple.h>
#include <io/blockinput.h>
#include <io/topology/in_topology.h>

#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <util/create_simulation.h>
#include <algorithm/create_md_sequence.h>

#include <time.h>

#include <config.h>

#include "check.h"
#include "check_forcefield.h"

using namespace std;


double finite_diff(topology::Topology & topo, 
		   configuration::Configuration &conf, 
		   simulation::Simulation & sim, 
		   interaction::Interaction &term,
		   size_t atom, size_t coord,
		   double const epsilon)
{
  conf.current().pos(atom)(coord) += epsilon;
  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim);
  conf.current().energies.calculate_totals();

  double e1 = conf.current().energies.potential_total;

  conf.current().pos(atom)(coord) -= 2 * epsilon;
  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim);
  conf.current().energies.calculate_totals();

  double e2 = conf.current().energies.potential_total;

  conf.current().pos(atom)(coord) += epsilon;

  return (e2 - e1) / 2.0 / epsilon;

}

int check_lambda_derivative(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    interaction::Interaction &term,
			    double const epsilon,
			    double const delta)
{
  int res, total = 0;
  
  CHECKING(term.name + " energy lambda derivative", res);

  // assume no different lambda dependence here!
  // this will be checked in a different function

  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives[0].zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives[0].calculate_totals();

  const double dE = conf.current().perturbed_energy_derivatives[0].potential_total;
  
  // change lambda, calculate energies
  topo.lambda(topo.lambda()+epsilon);
    
  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives[0].zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives[0].calculate_totals();

  double e1 = conf.current().energies.potential_total;

  // change lambda, calculate energies
  topo.lambda(topo.lambda()-2*epsilon);
 
  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives[0].zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives[0].calculate_totals();

  double e2 = conf.current().energies.potential_total;
  
  const double dfE = (e1 - e2) / (2 * epsilon);
  
  CHECK_APPROX_EQUAL(dE, dfE, delta, res);

  topo.lambda(topo.lambda() + epsilon);

  RESULT(res, total);

  return total;
  
  
}

int check_interaction(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      interaction::Interaction &term,
		      size_t atoms,
		      double const energy,
		      double const epsilon,
		      double const delta)
{
  int res, total = 0;
  
  CHECKING(term.name + " interaction energy", res);

  conf.current().force = 0;
  conf.current().energies.zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  CHECK_APPROX_EQUAL(conf.current().energies.potential_total,
		     energy, delta, res);
  RESULT(res, total);

  // finite diff
  CHECKING(term.name + " interaction force (finite diff)", res);

  for(size_t atom=0; atom < atoms; ++atom){

    conf.current().force = 0;
    conf.current().energies.zero();

    term.calculate_interactions(topo, conf, sim);
      
    math::Vec f = conf.current().force(atom);

    math::Vec finf;
    finf(0) = finite_diff(topo, conf, sim, term, atom, 0, epsilon);
    finf(1) = finite_diff(topo, conf, sim, term, atom, 1, epsilon);
    finf(2) = finite_diff(topo, conf, sim, term, atom, 2, epsilon);
	
    CHECK_APPROX_EQUAL(f(0), finf(0), delta, res);
    CHECK_APPROX_EQUAL(f(1), finf(1), delta, res);
    CHECK_APPROX_EQUAL(f(2), finf(2), delta, res);
  }
      
  RESULT(res, total);
  
  return total;

}


int check::check_forcefield(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    interaction::Forcefield & ff)
{
  int res=0, total=0;
  
  for(vector<interaction::Interaction *>::iterator
	it = ff.begin(),
	to = ff.end();
      it != to;
      ++it){
    
    // cout << (*it)->name << endl;

    if ((*it)->name == "QuarticBond"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 18.055276, 0.0000000001, 0.001);
    }
    else if ((*it)->name == "PerturbedQuarticBond"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 1.149568, 0.0000000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "Angle"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 12.170290, 0.00000001, 0.001);
    }
    else if ((*it)->name == "PerturbedAngle"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 0.714818, 0.00000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "ImproperDihedral"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 0.965060, 0.00000001, 0.001);      
    }
    else if ((*it)->name == "PerturbedImproperDihedral"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 2.642780, 0.00000001, 0.001);      
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "Dihedral"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 2.255206, 0.00000001, 0.001);    
    }
    else if ((*it)->name == "PerturbedDihedral"){
      total += check_interaction(topo, conf, sim, **it, topo.num_solute_atoms(), 13.314602, 0.00000001, 0.001);    
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "NonBonded"){
      total += check_interaction(topo, conf, sim, **it, topo.num_atoms(), -50.353066, 0.00000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else {
      CHECKING((*it)->name << " no check implemented!", res);
      CHECK_EQUAL(1, 2, res);
      RESULT(res, total);
    }
  }
  
  return total;
}
