/**
 * @file average.cc
 * implements Average member methods.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <math/volume.h>

#include "average.h"

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE energy

configuration::Average::Average()
  : Algorithm("Average"),
    m_sim_avg(),
    m_block_avg()
{
}

void configuration::Average::zero()
{
  m_sim_avg.zero();
  m_block_avg.zero();
}

int configuration::Average::apply(topology::Topology &topo, 
				  configuration::Configuration &conf,
				  simulation::Simulation &sim)
{

  DEBUG(7, "averaging...");
  
  m_sim_avg.update(conf.old().averages.simulation(), topo, conf, sim);

  if (sim.param().write.block_average)
    m_block_avg.update(conf.old().averages.block(), topo, conf, sim);
  
  return 0;

}

void configuration::Average::resize(topology::Topology const & topo,
				    configuration::Configuration const & conf,
				    simulation::Parameter const & param)
{
  DEBUG(10, "average: resize");

  m_sim_avg.resize(topo, conf, param);
  m_block_avg.resize(topo, conf, param);
  
}


void configuration::Average::Block_Average::zero()
{
  DEBUG(10, "Block_Average: zero");

  time = 0.0;

  energy_avg.zero();
  energy_fluct.zero();

  energy_derivative_avg.zero();
  energy_derivative_fluct.zero();
  
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      pressure_avg(a,b) = 0.0;
      pressure_fluct(a,b) = 0.0;
      virial_avg(a,b) = 0.0;
      virial_fluct(a,b) = 0.0;
      ekin_tensor_avg(a,b) = 0.0;
      ekin_tensor_fluct(a,b) = 0.0;

      box_avg(a)(b) = 0.0;
      box_fluct(a)(b) = 0.0;
    }
  }

  mass_avg = 0.0;
  mass_fluct = 0.0;
  
  temp_scaling_avg.assign(temp_scaling_avg.size(), 0.0);
  temp_scaling_fluct.assign(temp_scaling_avg.size(), 0.0);
  
  volume_avg = 0.0;
  volume_fluct = 0.0;

  lambda_avg = 0.0;
  lambda_fluct = 0.0;
  
}

void configuration::Average::Block_Average::
resize(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Parameter const & param)
{
  DEBUG(10, "Block_Average: resize");

  const unsigned int e_groups = unsigned(param.force.energy_group.size());
  const unsigned int baths = unsigned(param.multibath.multibath.size());

  energy_avg.resize(e_groups, baths);
  energy_fluct.resize(e_groups, baths);
  
  if (param.perturbation.perturbation){
    energy_derivative_avg.resize(e_groups, baths);
    energy_derivative_fluct.resize(e_groups, baths);
  }

  temp_scaling_avg.resize(baths);
  temp_scaling_fluct.resize(baths);

}

void configuration::Average::Block_Average::
update(Block_Average const & old,
      topology::Topology &topo, 
      configuration::Configuration &conf,
      simulation::Simulation &sim)
{

  DEBUG(7, "block averaging...");
  
  const double dt = sim.time_step_size();
  
  time = old.time + dt;

  update_energies(energy_avg,
		  energy_fluct,
		  conf.old().energies,
		  old.energy_avg,
		  old.energy_fluct,
		  dt,
		  0.0);
  
  if (sim.param().perturbation.perturbation){

    const double dlamt = sim.param().perturbation.dlamt;

    update_energies(energy_derivative_avg,
		    energy_derivative_fluct,
		    conf.old().perturbed_energy_derivatives,
		    old.energy_derivative_avg,
		    old.energy_derivative_fluct,
		    dt,
		    dlamt);
    
    lambda_avg = old.lambda_avg + topo.old_lambda() * dt;
    lambda_fluct = old.lambda_fluct + topo.old_lambda() * topo.old_lambda() * dt;

  }
  

  // and the rest...
  update_volumepressure(topo, conf, sim, old);

}

////////////////////////////////////////////////////
// accessors
////////////////////////////////////////////////////

void configuration::Average::Block_Average
::energy_average(configuration::Energy &energy, 
		 configuration::Energy &fluctuation)const
{
  energy_average_helper(energy_avg,
			energy_fluct,
			energy,
			fluctuation,
			0.0);
}

void configuration::Average::Block_Average
::energy_derivative_average(configuration::Energy &energy, 
			    configuration::Energy &fluctuation,
			    double & lambda,
			    double & lambda_fluct,
			    double const dlamt)const
{
  energy_average_helper(energy_derivative_avg,
			energy_derivative_fluct,
			energy,
			fluctuation,
			dlamt);
  
  lambda = lambda_avg / time;
  const double diff = this->lambda_fluct - lambda_avg * lambda_avg / time;
  if (diff > 0.0) lambda_fluct = sqrt(diff / time);
  else lambda_fluct = 0.0;
  
}

void configuration::Average::Block_Average
::pressure_average(math::Matrix &pressure, 
		   math::Matrix &pressure_fluctuations,
		   math::Matrix &virial,
		   math::Matrix &virial_fluctuations,
		   math::Matrix &kinetic_energy,
		   math::Matrix &kinetic_energy_fluctuations)const
{
  double diff;
  
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      pressure(a,b) = pressure_avg(a,b) / time;
      diff = pressure_fluct(a,b) -
	pressure_avg(a,b) * pressure_avg(a,b) / time;
      if (diff > 0.0)
	pressure_fluctuations(a,b) = sqrt(diff / time);
      else
	pressure_fluctuations(a,b) = 0.0;

      virial(a,b) = virial_avg(a,b) / time;
      diff = virial_fluct(a,b) -
	virial_avg(a,b) * virial_avg(a,b) / time;
      if (diff > 0.0)
	virial_fluctuations(a,b) = sqrt(diff / time);
      else
	virial_fluctuations(a,b) = 0.0;

      kinetic_energy(a,b) = ekin_tensor_avg(a,b) / time;
      diff = ekin_tensor_fluct(a,b) -
	ekin_tensor_avg(a,b) * ekin_tensor_avg(a,b) / time;
      if (diff > 0.0)
	kinetic_energy_fluctuations(a,b) = sqrt(diff / time);
      else
	kinetic_energy_fluctuations(a,b) = 0.0;

    }
  }

}

void configuration::Average::Block_Average
::mbs_average(double & mass, double & mass_fluctuations,
	      double & volume, double & volume_fluctuations,
	      math::Box & box, math::Box & box_fluctuations,
	      std::vector<double> & scaling,
	      std::vector<double> & scaling_fluctuations)const
{
  double diff;
  
  mass = mass_avg / time;
  diff = mass_fluct - mass_avg * mass_avg / time;
  if (diff > 0.0)
    mass_fluctuations = sqrt(diff / time);
  else
    mass_fluctuations = 0.0;

  volume = volume_avg / time;
  diff = volume_fluct - volume_avg * volume_avg / time;
  if (diff > 0.0)
    volume_fluctuations = sqrt(diff / time);
  else
    volume_fluctuations = 0.0;

  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      box(a)(b) = box_avg(a)(b) / time;
      diff = box_fluct(a)(b) -
	box_avg(a)(b) * box_avg(a)(b) / time;
      if (diff > 0.0)
	box_fluctuations(a)(b) = sqrt(diff / time);
      else
	box_fluctuations(a)(b) = 0.0;
    }
  }

  scaling.resize(temp_scaling_avg.size());
  scaling_fluctuations.resize(temp_scaling_avg.size());
  
  for(size_t i = 0; i < temp_scaling_avg.size(); ++i){

    scaling[i] = temp_scaling_avg[i] / time;
    diff = temp_scaling_fluct[i] -
      temp_scaling_avg[i] * temp_scaling_avg[i] / time;
    if (diff > 0.0)
      scaling_fluctuations[i] = sqrt(diff / time);
    else
      scaling_fluctuations[i] = 0.0;
  }
  
}



////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////

void configuration::Average::Block_Average
::update_energies(Energy &avg, Energy &fluct,
		  Energy const & e,
		  Energy const & old_avg,
		  Energy const & old_fluct,
		  double dt, double dlamt)
{

  double cumu = dt;
  if (dlamt) cumu = dlamt;
  
#define ENERGY_AVG(prop) \
avg.prop = old_avg.prop + \
cumu * e.prop; \
fluct.prop = old_fluct.prop + dt * e.prop * e.prop

  DEBUG(10, "average: old total = " << old_avg.total
	<< "\ttotal = " << e.total << "\tdt = " << dt << "\tdlamt = " << dlamt);

  // totals
  ENERGY_AVG(total);
  ENERGY_AVG(kinetic_total);
  ENERGY_AVG(potential_total);

  ENERGY_AVG(bonded_total);
  ENERGY_AVG(bond_total);
  ENERGY_AVG(angle_total);
  ENERGY_AVG(improper_total);
  ENERGY_AVG(dihedral_total);

  ENERGY_AVG(nonbonded_total);
  ENERGY_AVG(lj_total);
  ENERGY_AVG(crf_total);

  ENERGY_AVG(special_total);
  ENERGY_AVG(posrest_total);
  ENERGY_AVG(constraints_total);


  // kinetic energies...
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){

    assert(energy_avg.kinetic_energy.size() > i);
    assert(energy_fluct.kinetic_energy.size() > i);
    assert(e.kinetic_energy.size() > i);
    
    ENERGY_AVG(kinetic_energy[i]);

    assert(energy_avg.com_kinetic_energy.size() > i);
    assert(energy_fluct.com_kinetic_energy.size() > i);
    assert(e.com_kinetic_energy.size() > i);

    ENERGY_AVG(com_kinetic_energy[i]);

    assert(energy_avg.ir_kinetic_energy.size() > i);
    assert(energy_fluct.ir_kinetic_energy.size() > i);
    assert(e.ir_kinetic_energy.size() > i);

    ENERGY_AVG(ir_kinetic_energy[i]);

  }

  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){

    assert(energy_avg.bond_energy.size() > i);
    assert(energy_fluct.bond_energy.size() > i);

    ENERGY_AVG(bond_energy[i]);
    ENERGY_AVG(angle_energy[i]);
    ENERGY_AVG(improper_energy[i]);
    ENERGY_AVG(dihedral_energy[i]);

    // and the nonbonded groups
    assert(energy_avg.lj_energy.size() > i);
    assert(energy_fluct.lj_energy.size() > i);

    for(size_t j=0; j < e.lj_energy.size(); ++j){

      assert(energy_avg.lj_energy[i].size() > j);
      assert(energy_fluct.lj_energy[i].size() > j);

      ENERGY_AVG(lj_energy[i][j]);
      ENERGY_AVG(crf_energy[i][j]);
    }
  
    // special energies
    ENERGY_AVG(posrest_energy[i]);
    ENERGY_AVG(constraints_energy[i]);
  }

#undef ENERGY_AVG
}

void configuration::Average::Block_Average
::update_volumepressure(topology::Topology const & topo,
			configuration::Configuration const & conf,
			simulation::Simulation const & sim,
			Block_Average const & old)
{
  const double dt = sim.time_step_size();
  
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      // pressure
      pressure_avg(a,b) = old.pressure_avg(a,b) + 
	dt * conf.old().pressure_tensor(a,b);
      pressure_fluct(a,b) = old.pressure_fluct(a,b) +
	dt * conf.old().pressure_tensor(a,b) * conf.old().pressure_tensor(a,b);
      // virial
      virial_avg(a,b) = old.virial_avg(a,b) + 
	dt * conf.old().virial_tensor(a,b);
      virial_fluct(a,b) = old.virial_fluct(a,b) +
	dt * conf.old().virial_tensor(a,b) * conf.old().virial_tensor(a,b);
      // ekin_tensor
      ekin_tensor_avg(a,b) = old.ekin_tensor_avg(a,b) + 
	dt * conf.old().kinetic_energy_tensor(a,b);
      ekin_tensor_fluct(a,b) = old.ekin_tensor_fluct(a,b) +
	dt * conf.old().kinetic_energy_tensor(a,b) * 
	conf.old().kinetic_energy_tensor(a,b);

      // box
      box_avg(a)(b) = old.box_avg(a)(b) + 
	dt * conf.old().box(a)(b);
      box_fluct(a)(b) = old.box_fluct(a)(b) +
	dt * conf.old().box(a)(b) * 
	conf.old().box(a)(b);

    }
  }

  const double tot_mass = math::sum(topo.mass());
  mass_avg = old.mass_avg + dt * tot_mass;
  mass_fluct = old.mass_fluct + dt * tot_mass * tot_mass;
  const double volume = math::volume(conf.old().box, conf.boundary_type);
  volume_avg = old.volume_avg + dt * volume;
  volume_fluct = old.volume_fluct + dt * volume * volume;
  
  simulation::Multibath const & m = sim.multibath();

  for(size_t i=0; i < temp_scaling_avg.size(); ++i){

    temp_scaling_avg[i] = old.temp_scaling_avg[i] + 
      dt * m[i].scale;
    temp_scaling_fluct[i] = old.temp_scaling_fluct[i] +
      dt * m[i].scale * m[i].scale;

  }

}

void configuration::Average::Block_Average
::energy_average_helper(configuration::Energy const & avg,
			configuration::Energy const & fluct,
			configuration::Energy &energy, 
			configuration::Energy &fluctuation,
			double const dlamt)const
{
  DEBUG(10, "new average total: " << energy_avg.total << "\ttime: " << time);

  double diff;
  double cumu = time;
  if (dlamt) cumu = 1;
  
  Energy &e = energy;
  Energy &f = fluctuation;

  e.resize(unsigned(energy_avg.bond_energy.size()),
		unsigned(energy_avg.kinetic_energy.size()));
  f.resize(unsigned(energy_avg.bond_energy.size()),
	  unsigned(energy_avg.kinetic_energy.size()));

#define ENERGY_RES(prop) \
  e.prop = avg.prop / cumu; \
  diff = fluct.prop - avg.prop * avg.prop / time; \
  if (diff > 0.0) f.prop = sqrt(diff / time); \
  else f.prop = 0.0
  
  ENERGY_RES(total);
  ENERGY_RES(kinetic_total);
  ENERGY_RES(potential_total);

  ENERGY_RES(bonded_total);
  ENERGY_RES(bond_total);
  ENERGY_RES(angle_total);
  ENERGY_RES(improper_total);
  ENERGY_RES(dihedral_total);

  ENERGY_RES(nonbonded_total);
  ENERGY_RES(lj_total);
  ENERGY_RES(crf_total);

  ENERGY_RES(special_total);
  ENERGY_RES(posrest_total);
  ENERGY_RES(constraints_total);

  // std::cout << e.kinetic_energy.size() << std::endl;

  // kinetic energies...
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){

    ENERGY_RES(kinetic_energy[i]);
    ENERGY_RES(com_kinetic_energy[i]);
    ENERGY_RES(ir_kinetic_energy[i]);

  }

  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){

    assert(energy_avg.bond_energy.size() > i);
    assert(energy_fluct.bond_energy.size() > i);

    ENERGY_RES(bond_energy[i]);
    ENERGY_RES(angle_energy[i]);
    ENERGY_RES(improper_energy[i]);
    ENERGY_RES(dihedral_energy[i]);

    // and the nonbonded groups
    assert(energy_avg.lj_energy.size() > i);
    assert(energy_fluct.lj_energy.size() > i);

    for(size_t j=0; j < e.lj_energy.size(); ++j){

      assert(energy_avg.lj_energy[i].size() > j);
      assert(energy_fluct.lj_energy[i].size() > j);

      ENERGY_RES(lj_energy[i][j]);
      ENERGY_RES(crf_energy[i][j]);
    }
  
    // special energies
    ENERGY_RES(posrest_energy[i]);
    ENERGY_RES(constraints_energy[i]);
  }

#undef ENERGY_RES
}
