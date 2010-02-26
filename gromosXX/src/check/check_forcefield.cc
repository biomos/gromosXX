/**
 * @file check_forcefield.cc
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/nonbonded/interaction/nonbonded_interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/standard_pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/extended_grid_pairlist_algorithm.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

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
		   double const epsilon,
		   bool physical=true)
{
  conf.current().pos(atom)(coord) += epsilon; //epsilon from end of file check::check_forcefield
  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim); // from interaction/interaction.h set to 0?
  conf.current().energies.calculate_totals(); // from configuration/energy.cc

  double e1;
  if(physical)
    e1= conf.current().energies.potential_total;
  else
    e1= conf.current().energies.special_total;

  conf.current().pos(atom)(coord) -= 2 * epsilon;
  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim);
  conf.current().energies.calculate_totals();

  double e2;
  if(physical)
    e2= conf.current().energies.potential_total;
  else
    e2= conf.current().energies.special_total;

  conf.current().pos(atom)(coord) += epsilon;

  // std::cout << "atom=" << atom << " e1=" << e1 
  // << " e2=" << e2 << " epsilon=" << epsilon << std::endl;

  return (e2 - e1) / 2.0 / epsilon;

}

int nonbonded_hessian(topology::Topology & topo, 
		      configuration::Configuration &conf, 
		      simulation::Simulation & sim, 
		      interaction::Nonbonded_Interaction &term,
		      size_t atom_i, size_t atom_j,
		      double const epsilon,
		      double const delta,
		      int & res)
{
//hessian matrix is square matrix of second-order partial derivatives of a function
  math::Vec pos_i = conf.current().pos(atom_i);
  math::Vec f1, f2; //from algorithm/constraints/flexible_constraint.cc they are hessian functions
  double e_lj, e_crf;
  math::Matrix fd_h, h;

  for(int d=0; d<3; ++d){
    
    conf.current().pos(atom_i)(d) += epsilon;

    term.calculate_interaction(topo, conf, sim, atom_i, atom_j, f1, e_lj, e_crf); //calculate interaction for a given atom pair

    conf.current().pos(atom_i)(d) -= 2 * epsilon;
    
    term.calculate_interaction(topo, conf, sim, atom_i, atom_j, f2, e_lj, e_crf);

    // restore...
    conf.current().pos(atom_i) = pos_i;

    for(int i=0; i<3; ++i)
      fd_h(i, d) = (f2(i) - f1(i)) / 2.0 / epsilon;

  }

  conf.current().force = 0;
  conf.current().energies.zero();

  // need to create a pairlist
  term.calculate_interactions(topo, conf, sim);
  
  term.calculate_hessian(topo, conf, sim, atom_i, atom_j, h);

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      CHECK_APPROX_EQUAL(h(i,j), fd_h(i,j), delta, res);

  return res;

}


int check_lambda_derivative(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    interaction::Interaction &term,
			    double const epsilon,
			    double const delta, //from end of file check::check_forcefield
			    bool physical=true
)
{
  int res, total = 0; // checking how many test, how many negative results
  
  std::string name = term.name;

  if (term.name == "NonBonded"){

    if (sim.param().force.interaction_function == simulation::cgrain_func)
      name = "CG-" + name;    
    
    if (sim.param().force.special_loop == simulation::special_loop_spc)
      if (sim.param().pairlist.grid)
	name += " (grid, spc loop)";
      else
	name += " (std, spc loop)";
    else
      if (sim.param().pairlist.grid)
	name += " (grid, std loop)";
      else
	name += " (std, std loop)";
  }

  CHECKING(name + " energy lambda derivative", res);

  // assume no different lambda dependence here!
  // this will be checked in a different function

  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives.zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives.calculate_totals();

  double dE;
  if (physical) dE = conf.current().perturbed_energy_derivatives.potential_total;
  else dE = conf.current().perturbed_energy_derivatives.special_total;
  
  // change lambda, calculate energies
  topo.lambda(topo.lambda()+epsilon);
  topo.update_for_lambda();
  
  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives.zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives.calculate_totals();

  double e1;
  if(physical)
    e1 = conf.current().energies.potential_total;
  else 
    e1 = conf.current().energies.special_total;

  // change lambda, calculate energies
  topo.lambda(topo.lambda()-2*epsilon);
  topo.update_for_lambda();
  
  conf.current().force = 0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives.zero();

  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  conf.current().perturbed_energy_derivatives.calculate_totals();

  double e2;
  if(physical)
    e2 = conf.current().energies.potential_total;
  else
    e2 = conf.current().energies.special_total;


  const double dfE = (e1 - e2) / (2 * epsilon);
  
  CHECK_APPROX_EQUAL(dE, dfE, delta, res);

  // std::cout << " dE " << dE << " dfE " << dfE << std::endl;
  
  topo.lambda(topo.lambda() + epsilon);
  topo.update_for_lambda();
  
  RESULT(res, total);

  return total;
  
  
}

int check_interaction(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      interaction::Interaction &term,
		      size_t atoms,
		      double const energy, //from end of file check::check_forcefield
		      double const epsilon,
		      double const delta)
{
  int res, total = 0;

  std::string name = term.name;

  if (term.name == "NonBonded"){
    if (sim.param().force.interaction_function == simulation::cgrain_func)
      name = "CG-" + name;
    if (sim.param().nonbonded.rf_excluded == 1)
      name = "newRF-" + name ;
    if (sim.param().force.special_loop == simulation::special_loop_spc)
      if (sim.param().pairlist.grid)
	name += " (grid, spc loop)";
      else
	name += " (std, spc loop)";
    else
      if (sim.param().pairlist.grid)
	name += " (grid, std loop)";
      else
	name += " (std, std loop)";
  }

  CHECKING(name + " interaction energy", res);
  
  conf.current().force = 0;
  conf.current().energies.zero();

  term.calculate_interactions(topo, conf, sim);
  
  conf.current().energies.calculate_totals();
  CHECK_APPROX_EQUAL(conf.current().energies.potential_total,
		     energy, delta, res);
  RESULT(res, total);

  // finite diff
  CHECKING(name + " interaction force (finite diff)", res);
   
  for(size_t atom=0; atom < atoms; ++atom){
 
    conf.current().force = 0;
    conf.current().energies.zero();

    term.calculate_interactions(topo, conf, sim);
   
    math::Vec f = conf.current().force(atom);
 //   cout << f(0) << " " << f(1) << " " << f(2) << endl;
    math::Vec finf;
    finf(0) = finite_diff(topo, conf, sim, term, atom, 0, epsilon,true); // 0 is coord x,y or z
    finf(1) = finite_diff(topo, conf, sim, term, atom, 1, epsilon,true);
    finf(2) = finite_diff(topo, conf, sim, term, atom, 2, epsilon,true);
  //  cout << finf(0) << " " << finf(1) << " " << finf(2) << endl;
   
    CHECK_APPROX_EQUAL(f(0), finf(0), delta, res);
    CHECK_APPROX_EQUAL(f(1), finf(1), delta, res);
    CHECK_APPROX_EQUAL(f(2), finf(2), delta, res);
  }
      
  RESULT(res, total);
  
  if (term.name == "NonBonded" && 
      sim.param().force.interaction_function != simulation::cgrain_func && 
      sim.param().nonbonded.rf_excluded == 0){
    
    // std::cout << "num atoms : " << topo.num_atoms() << std::endl;
    // std::cout << "num solvents : " << topo.num_solvents() << std::endl;
    // std::cout << "num solvent atoms : " << topo.num_solvent_atoms() << std::endl;

    interaction::Nonbonded_Interaction & ni =
      dynamic_cast<interaction::Nonbonded_Interaction &>(term);

    CHECKING(name + " hessian (finite diff)", res);
    res += nonbonded_hessian(topo, conf, sim, ni, 1, 8, epsilon, delta, res);
    res += nonbonded_hessian(topo, conf, sim, ni, 2, 11, epsilon, delta, res);
    res += nonbonded_hessian(topo, conf, sim, ni, 9, 69, epsilon, delta, res);
    res += nonbonded_hessian(topo, conf, sim, ni, 9, 70, epsilon, delta, res);
    RESULT(res, total);
  }

  return total;

}

int check_distanceres_interaction(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       interaction::Interaction &term,
			       size_t atoms,
			       double const energy,
			       double const epsilon,
			       double const delta)
{
  int res, total = 0;
  
  std::string name = term.name;
  
  CHECKING(name + " interaction energy", res);

  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  CHECK_APPROX_EQUAL(conf.current().energies.special_total,
		     energy, delta, res);
  RESULT(res, total);
  
  // finite diff
  CHECKING(name + " interaction force (finite diff)", res);
  
  for(size_t i=0; i < topo.distance_restraints().size(); ++i){

    int atom = topo.distance_restraints()[i].v1.atom(0);
    
    conf.current().force = 0;
    conf.current().energies.zero();

    term.calculate_interactions(topo, conf, sim);
      
    math::Vec f = conf.current().force(atom);

    math::Vec finf;
    finf(0) = finite_diff(topo, conf, sim, term, atom, 0, epsilon, false); 
    finf(1) = finite_diff(topo, conf, sim, term, atom, 1, epsilon, false);
    finf(2) = finite_diff(topo, conf, sim, term, atom, 2, epsilon, false);
	//finite_diff from upwards
    CHECK_APPROX_EQUAL(f(0), finf(0), delta, res);
    CHECK_APPROX_EQUAL(f(1), finf(1), delta, res);
    CHECK_APPROX_EQUAL(f(2), finf(2), delta, res);
  }
      
  RESULT(res, total);

  return total;

}
int check_dihrest_interaction(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      interaction::Interaction &term,
			      size_t atoms,
			      double const energy,
			      double const epsilon,
			      double const delta)
{
  int res, total = 0;
  
  std::string name = term.name;
  
  CHECKING(name + " interaction energy", res);

  conf.current().force = 0;
  conf.current().energies.zero();
  
  term.calculate_interactions(topo, conf, sim);
      
  conf.current().energies.calculate_totals();
  CHECK_APPROX_EQUAL(conf.current().energies.special_total,
		     energy, delta, res);
  RESULT(res, total);
  
  // finite diff
  CHECKING(name + " interaction force (finite diff)", res);
  
  for(size_t i=0; i < topo.dihedral_restraints().size(); ++i){

    std::vector<int> atom;
    atom.push_back(topo.dihedral_restraints()[i].i);
    atom.push_back(topo.dihedral_restraints()[i].j);
    atom.push_back(topo.dihedral_restraints()[i].k);
    atom.push_back(topo.dihedral_restraints()[i].l);
    
    conf.current().force = 0;
    conf.current().energies.zero();

    term.calculate_interactions(topo, conf, sim);
    
    for(size_t j=0; j<atom.size(); ++j){
	
      math::Vec f = conf.current().force(atom[j]);
      
      math::Vec finf;
      
      finf(0) = finite_diff(topo, conf, sim, term, atom[j], 0, epsilon, false);
      finf(1) = finite_diff(topo, conf, sim, term, atom[j], 1, epsilon, false);
      finf(2) = finite_diff(topo, conf, sim, term, atom[j], 2, epsilon, false);
      
      CHECK_APPROX_EQUAL(f(0), finf(0), delta, res);
      CHECK_APPROX_EQUAL(f(1), finf(1), delta, res);
      CHECK_APPROX_EQUAL(f(2), finf(2), delta, res);
    }
  }
  
      
  RESULT(res, total);

  return total;

}


int check::check_forcefield(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    interaction::Forcefield & ff,
			    std::map<std::string, double> & ref)
{
  int res=0, total=0;
  
  for(vector<interaction::Interaction *>::iterator
	it = ff.begin(),  //pointer to begin of a c++ vector
	to = ff.end();
      it != to;
      ++it){
    
    // cout << (*it)->name << endl;

    if ((*it)->name == "QuarticBond"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(),
				 ref["QuarticBond"], 
				 0.0000000001, 0.001);
    }
    else if ((*it)->name == "PerturbedQuarticBond"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["PerturbedQuarticBond"], 
				 0.0000000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "Angle"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["Angle"],
				 0.00000001, 0.001);
    }
    else if ((*it)->name == "PerturbedAngle"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["PerturbedAngle"],
				 0.00000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "ImproperDihedral"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["ImproperDihedral"],
				 0.00000001, 0.001);      
    }
    else if ((*it)->name == "PerturbedImproperDihedral"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["PerturbedImproperDihedral"],
				 0.00000001, 0.001);      
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "Dihedral"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["Dihedral"],
				 0.00000001, 0.001);    
    }
    else if ((*it)->name == "Crossdihedral"){
      total += check_interaction(topo, conf, sim, **it,
				 topo.num_solute_atoms(),
				 ref["Crossdihedral"],
				 0.00000001, 0.001);
    }
    else if ((*it)->name == "PerturbedDihedral"){
      total += check_interaction(topo, conf, sim, **it, 
				 topo.num_solute_atoms(), 
				 ref["PerturbedDihedral"],
				 0.00000001, 0.001);    
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
    }
    else if ((*it)->name == "NonBonded"){
      if (sim.param().force.interaction_function == simulation::cgrain_func){

	sim.param().force.special_loop = simulation::special_loop_off;
	total += check_interaction(topo, conf, sim, **it, topo.num_atoms(), 
				   ref["NonBonded_cg"],
				   0.0000000001, 0.01);
	total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);

      }
      else{
	sim.param().force.special_loop = simulation::special_loop_spc;
	total += check_interaction(topo, conf, sim, **it, topo.num_atoms(), 
				   ref["NonBonded"], 0.00000001, 0.001);
	total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
	sim.param().force.special_loop = simulation::special_loop_off;
	//total += check_interaction(topo, conf, sim, **it, topo.num_atoms(), 
	//			   ref["NonBonded"], 0.00000001, 0.001);
	//total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001);
	if (! sim.param().pairlist.grid){ //if not 1, eg grid see in_parameter.cc
	  // construct a "standard" pairlist algorithm 
	  interaction::Pairlist_Algorithm * pa = new interaction::Standard_Pairlist_Algorithm;
	  interaction::Nonbonded_Interaction * ni = dynamic_cast<interaction::Nonbonded_Interaction *>(*it);
	  
	  pa->set_parameter(&ni->parameter());
	  interaction::Pairlist_Algorithm * old_pa = &ni->pairlist_algorithm();
	  ni->pairlist_algorithm(pa);
	  
	  sim.param().pairlist.grid = false;
	  
	  ni->init(topo, conf, sim, std::cout, true);
          
	  total += check_interaction(topo, conf, sim, *ni, topo.num_atoms(), 
				     ref["NonBonded"], 0.00000001, 0.001);
	  total += check_lambda_derivative(topo, conf, sim, *ni, 0.001, 0.001);
	  
        //set RF calculation for excluded pairs on
          sim.param().nonbonded.rf_excluded = 1;
          total += check_interaction(topo, conf, sim, **it, topo.num_atoms(), 
				   ref["NonBonded_newRF"], 0.00000001, 0.001);
	  sim.param().pairlist.grid = true;
          sim.param().nonbonded.rf_excluded = 0;
          delete pa;
          // construct a "grid" pairlist algorithm 
       //   std::cout << "before new" << std::endl;
	  pa = new interaction::Extended_Grid_Pairlist_Algorithm;
       //   std::cout << "nre" << std::endl;
	  pa->set_parameter(&ni->parameter());
          sim.param().pairlist.grid = true;
	  sim.param().pairlist.grid_size = 0.1;
         // std::cout << "before init" << std::endl;
          ni->pairlist_algorithm(pa);
	  ni->init(topo, conf, sim, std::cout, true);
	//  std::cout << "before check" << std::endl;
	  total += check_interaction(topo, conf, sim, *ni, topo.num_atoms(), 
				     ref["NonBonded"], 0.00000001, 0.001);
	  total += check_lambda_derivative(topo, conf, sim, *ni, 0.001, 0.001);
	//  std::cout << "after check" << std::endl;
	  sim.param().pairlist.grid = false;
	  ni->pairlist_algorithm(old_pa);
	  ni->init(topo, conf, sim, std::cout, true);
	  delete pa;
          
	}
	else{
	  // construct a "grid" pairlist algorithm 
	  interaction::Pairlist_Algorithm * pa = new interaction::Extended_Grid_Pairlist_Algorithm;
	  interaction::Nonbonded_Interaction * ni = dynamic_cast<interaction::Nonbonded_Interaction *>(*it);
	  
	  pa->set_parameter(&ni->parameter());
	  interaction::Pairlist_Algorithm * old_pa = &ni->pairlist_algorithm();
	  ni->pairlist_algorithm(pa);
	  
	  sim.param().pairlist.grid = true;
	  sim.param().pairlist.grid_size = 0.1;
	  ni->init(topo, conf, sim, std::cout, true);
	  
	  total += check_interaction(topo, conf, sim, *ni, topo.num_atoms(), 
				     ref["NonBonded"], 0.00000001, 0.001);
	  total += check_lambda_derivative(topo, conf, sim, *ni, 0.001, 0.001);
	  
	  sim.param().pairlist.grid = false;
	  ni->pairlist_algorithm(old_pa);
	  ni->init(topo, conf, sim, std::cout, true);
	  delete pa;
	}
      }
    }

    else if ((*it)->name == "DistanceRestraint"){
      total += check_distanceres_interaction(topo, conf, sim, **it, 
					     topo.num_solute_atoms(), 
					     ref["DistanceRestraint"],
					     0.0000000001, 0.001);
    }
    else if ((*it)->name == "PerturbedDistanceRestraint"){
      total += check_distanceres_interaction(topo, conf, sim, **it, 
					     topo.num_solute_atoms(), 
					     ref["PerturbedDistanceRestraint"],
					     0.0000000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001, false);
    }
    else if ((*it)->name == "DihedralRestraint"){
      total += check_dihrest_interaction(topo, conf, sim, **it, 
					 topo.num_solute_atoms(), 
					 ref["DihedralRestraint"],
					 0.0000000001, 0.001);
    }
    else if ((*it)->name == "PerturbedDihedralRestraint"){
      total += check_dihrest_interaction(topo, conf, sim, **it, 
					 topo.num_solute_atoms(), 
					 ref["PerturbedDihedralRestraint"],
					 0.0000000001, 0.001);
      total += check_lambda_derivative(topo, conf, sim, **it, 0.001, 0.001, false);
    }
    else if ((*it)->name == "MolecularVirial"){
      // no real interaction...
    }
    else {
      CHECKING((*it)->name << " no check implemented!", res);
      CHECK_EQUAL(1, 2, res);
      RESULT(res, total);
    }
  }
  
  return total;
}

int check::check_atomic_cutoff(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       interaction::Forcefield & ff,
			       std::map<std::string, double> & ref)
{
  int res=0, total=0;
  
  for(vector<interaction::Interaction *>::iterator
	it = ff.begin(),
	to = ff.end();
      it != to;
      ++it){
    
    if ((*it)->name == "NonBonded"){
  
      CHECKING((*it)->name + ": atomic cutoff", res);

      sim.param().force.special_loop = simulation::special_loop_off;
      
      conf.current().force = 0;
      conf.current().energies.zero();

      (*it)->calculate_interactions(topo, conf, sim);
      
      conf.current().energies.calculate_totals();
      CHECK_APPROX_EQUAL(conf.current().energies.potential_total,
		         ref["NonBonded_atomic"], 0.001, res);
      RESULT(res, total);
    }
  }
  
  return total;
}
