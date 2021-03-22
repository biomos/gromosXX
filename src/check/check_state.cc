/**
 * @file check_forcefield.cc
 */


#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/error.h"

#include "../interaction/interaction_types.h"
#include "../io/instream.h"
#include "../util/parse_tcouple.h"
#include "../io/blockinput.h"
#include "../io/topology/in_topology.h"

#include "../algorithm/integration/leap_frog.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/berendsen_thermostat.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"

#include "../interaction/forcefield/create_forcefield.h"

#include "../util/create_simulation.h"
#include "../algorithm/create_md_sequence.h"

#include "../math/periodicity.h"
#include "../math/volume.h"

#include "../util/template_split.h"
#include "../util/prepare_virial.h"

#include <time.h>

#include "check.h"
#include "check_state.h"

using namespace std;

template<math::boundary_enum b>
void scale_positions(topology::Topology & topo,
		     configuration::Configuration & conf,
                     simulation::Simulation & sim,
		     math::Vec const scale)
{
  // calculate rel_mol_com_pos ;-)
  math::Periodicity<b> periodicity(conf.current().box);

  std::vector<math::Vec> com_pos; //center of mass position
  std::vector<math::Matrix> com_ekin;
  
  util::centre_of_mass(topo, conf, com_pos, com_ekin, sim);

  topology::Pressuregroup_Iterator
    pg_it = topo.pressure_group_begin(),
    pg_to = topo.pressure_group_end();

  for(int i=0; pg_it != pg_to; ++pg_it, ++i){
    
    topology::Atom_Iterator a_it = pg_it.begin(),
      a_to = pg_it.end();

    math::VArray &pos = conf.current().pos;

    for( ; a_it != a_to; ++a_it){
      math::Vec v;
      periodicity.nearest_image(pos(*a_it), com_pos[i], v);

      for(int i=0; i<3; ++i)
	pos(*a_it)(i) = (pos(*a_it)(i) - v(i)) * scale(i) + v(i);
    }
  }

  for(unsigned int i=0; i<3; ++i){
    conf.current().box(0)(i) *=scale(i);
    conf.current().box(1)(i) *=scale(i);
    conf.current().box(2)(i) *=scale(i);
  }
  
}

void scale_positions_atomic(topology::Topology & topo,
			    configuration::Configuration & conf,
			    math::Vec const scale)
{
  for(size_t i = 0; i< topo.num_atoms(); ++i){
    for(unsigned int j=0; j<3; ++j)
      conf.current().pos(i)(j) *= scale(j);
  }
  
  for(unsigned int i=0; i<3; ++i){
    conf.current().box(0)(i) *=scale(i);
    conf.current().box(1)(i) *=scale(i);
    conf.current().box(2)(i) *=scale(i);
  }
}


//================================================================================
// GENERAL STATE CHECKING (molecular virial)
//================================================================================

int check::check_state(topology::Topology & mytopo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       interaction::Forcefield & ff)
{
  const double epsilon = 0.000001;
  int res=0, total=0;

  CHECKING("molecular virial (finite diff)",res);

  math::Box box_backup = conf.current().box;
  math::VArray pos_backup = conf.current().pos;
  topology::Topology topo(mytopo);

  conf.exchange_state();
  conf.current().pos = conf.old().pos;
  conf.current().box = conf.old().box;

  for(size_t s = 0; s < topo.num_solute_atoms(); ++s){
    topo.one_four_pair(s).clear();
    for(size_t t=s+1; t < topo.num_solute_atoms(); ++t){
      
      topo.all_exclusion(s).insert(t);
    }
  }

  // prepare rel_mol_com_pos
  util::prepare_virial(topo, conf, sim);

  // calculate the virial (real)
  ff.apply(topo, conf, sim);

  conf.exchange_state();
  
  algorithm::Pressure_Calculation * pcalc =
    new algorithm::Pressure_Calculation;

  pcalc->apply(topo, conf, sim);

  delete pcalc;

  // The virial is internally stored without the factor -0.5
  math::Matrix realP(-0.5*conf.old().virial_tensor);

  // std::cout << "real virial is\n"
  // << math::m2s(realP)
  // << "\n--------------------------------------------------\n" << std::endl;
  
  math::Matrix finP;
  finP=0;

  double rmsf = 0.0;
  if (sim.param().force.interaction_function == simulation::lj_ls_func)
    rmsf = sim.param().nonbonded.influence_function_rms_force_error;
  
  for(int i=0; i < 3; i++){
    math::Vec s1(1); 
    s1(i) += epsilon;
      
    SPLIT_BOUNDARY(scale_positions, topo, conf, sim, s1); //assign boundaries
    ff.apply(topo, conf, sim);
    conf.current().energies.calculate_totals();
    
    double e1 = conf.current().energies.nonbonded_total;
      
    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;

    math::Vec s2(1);
    s2(i) -= epsilon; 

    SPLIT_BOUNDARY(scale_positions, topo, conf, sim, s2);
    ff.apply(topo, conf, sim);
    conf.current().energies.calculate_totals();

    double e2 = conf.current().energies.nonbonded_total;

    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;

    finP(i,i) = -0.5 * conf.current().box(i)(i) * (e2-e1)/(2*epsilon) /
      conf.current().box(i)(i);
  }

  for(int i=0; i < 3; i++){

    CHECK_APPROX_EQUAL_RMSFERR(realP(i,i), finP(i,i), 0.00001, rmsf, res);
    if (res) {
    std::cout << "is: " << std::setw(15) << realP(i,i)
            << " should: " << setw(15) << finP(i,i) << std::endl;
    }
  }

  RESULT(res,total);

  conf.current().box = box_backup;
  conf.current().pos = pos_backup;

  return total;
}

//================================================================================
// ATOMIC VIRIAL CHECKING
//================================================================================

static void zero(configuration::Configuration & conf)
{
  conf.current().force = 0.0;
  conf.current().energies.zero();
  conf.current().perturbed_energy_derivatives.zero();

  conf.current().virial_tensor = 0.0;
}

static void fdiff_atomic_virial(topology::Topology & topo, 
				configuration::Configuration &conf, 
				simulation::Simulation & sim, 
				interaction::Interaction &term,
				math::Vec &vir,
				double const epsilon)
{
  for(int i=0; i<3; ++i){

    math::Vec scale(1);
    scale(i) += epsilon;
    scale_positions_atomic(topo, conf, scale);

    zero(conf);
    term.calculate_interactions(topo, conf, sim);
    conf.current().energies.calculate_totals();
    const double e1 = conf.current().energies.potential_total;
    
    // revert positions
    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;
    
    // do the fun backwards
    scale = 1;
    scale(i) -= epsilon;
    scale_positions_atomic(topo, conf, scale);

    zero(conf);
    term.calculate_interactions(topo, conf, sim);
    conf.current().energies.calculate_totals();
    const double e2 = conf.current().energies.potential_total;
    
    // revert
    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;
    
    vir(i) = -0.5 * conf.current().box(i)(i) * (e2-e1)/(2*epsilon) /
      conf.current().box(i)(i);
  }
  
}

static void exact_atomic_virial(topology::Topology & topo, 
				configuration::Configuration &conf, 
				simulation::Simulation & sim, 
				interaction::Interaction &term,
				algorithm::Pressure_Calculation *pcalc,
				math::Vec &vir,
				double const accuracy,
				int & res,
				int & total)
{
  zero(conf);
  term.calculate_interactions(topo, conf, sim);
  // get energies, virial into the correct state (w/o leap frog)
  // positions are not affected
  conf.exchange_state();

  // calculate the pressure...
  pcalc->apply(topo, conf, sim);

  conf.exchange_state();
 
  // the virial is internally stored without the prefactor of -0.5  
  math::Matrix realP(-0.5*conf.current().virial_tensor);

  for(int i=0; i < 3; i++){
    // check the trace
    CHECK_APPROX_EQUAL(realP(i,i), vir(i), accuracy, res);
  }
  
  RESULT(res,total);
  
}

int check::check_atomic_virial(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       interaction::Forcefield & ff)
{

  //set RF calculation for excluded pairs off 
  //(became necessary after removal of GROMOS96COMPAT block)
  sim.param().nonbonded.rf_excluded = 0; 

  const double epsilon = 0.0000001;
  const double accuracy = 0.00001;
  
  int res=0, total=0;

  math::Vec vir;

  algorithm::Pressure_Calculation * pcalc =
    new algorithm::Pressure_Calculation;

  conf.exchange_state();
  conf.current().pos = conf.old().pos;
  conf.current().box = conf.old().box;

  // pick out the correct interaction
  for(vector<interaction::Interaction *>::iterator
	it = ff.begin(),
	to = ff.end();
      it != to;
      ++it){
    
    CHECKING("atomic virial ("+(*it)->name +")",res);
    
    fdiff_atomic_virial(topo, conf, sim, **it, vir, epsilon);
    exact_atomic_virial(topo, conf, sim, **it, pcalc, vir, accuracy, res, total);
    
  } // loop over the forcefield terms

  return total;
}
