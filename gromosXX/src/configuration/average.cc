/**
 * @file average.cc
 * implements Average member methods.
 */

#include "../stdheader.h"
#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../math/volume.h"

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

  // set all sasa and volume values to zero
  sasatot_avg = 0.0;
  sasatot_fluct = 0.0;

  sasa_buriedvol_tot_avg = 0.0;
  sasa_buriedvol_tot_fluct = 0.0;
  
  sasa_avg.assign(sasa_avg.size(), 0.0);
  sasa_fluct.assign(sasa_avg.size(), 0.0);
  
  sasa_buriedvol_avg.assign(sasa_buriedvol_avg.size(), 0.0);
  sasa_buriedvol_fluct.assign(sasa_buriedvol_avg.size(), 0.0);
  
}

void configuration::Average::Block_Average::
resize(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Parameter const & param)
{
  DEBUG(10, "Block_Average: resize");

  const unsigned int e_groups = unsigned(param.force.energy_group.size());
  const unsigned int baths = unsigned(param.multibath.multibath.size());
  const unsigned int num_sasa_atoms = unsigned(topo.sasa_parameter().size());
// ANITA
  const unsigned int nr_lambdas = unsigned(param.precalclam.nr_lambdas);

//  energy_avg.resize(e_groups, baths);
//  energy_fluct.resize(e_groups, baths);
  energy_avg.resize(e_groups, baths, nr_lambdas);
  energy_fluct.resize(e_groups, baths, nr_lambdas);
  
  if (param.perturbation.perturbation){
//    energy_derivative_avg.resize(e_groups, baths);
//    energy_derivative_fluct.resize(e_groups, baths);
    energy_derivative_avg.resize(e_groups, baths, nr_lambdas);
    energy_derivative_fluct.resize(e_groups, baths,nr_lambdas);
  }

  temp_scaling_avg.resize(baths);
  temp_scaling_fluct.resize(baths);

  // resize sasa and volume vectors
  if (param.sasa.switch_sasa) {
    sasa_avg.resize(num_sasa_atoms);
    sasa_fluct.resize(num_sasa_atoms);

    sasa_buriedvol_avg.resize(num_sasa_atoms);
    sasa_buriedvol_fluct.resize(num_sasa_atoms);
  }
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

  DEBUG(5, "ANITA, Block_Average::update, checking nr lambdas" 
            << "\n conf.old: " << conf.old().energies.A_lj_total.size());

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

  // update sasa and volume averages/fluctuations
  if (sim.param().sasa.switch_sasa) {
    update_sasa(topo, conf, sim, old);
    if (sim.param().sasa.switch_volume)
      update_sasavol(topo, conf, sim, old);
  }

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
  double diff = 0.0;
  
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
  double diff = 0.0;
  
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

void configuration::Average::Block_Average
::sasa_average(std::vector<double> & sasa,
               std::vector<double> & sasa_fluctuations,
               double & sasatot, double & sasatot_fluct)const
{
  sasa.resize(sasa_avg.size());
  sasa_fluctuations.resize(sasa_avg.size());
  
  for (unsigned int i = 0; i < sasa_avg.size(); ++i){
    sasa[i] = sasa_avg[i] / time;
    const double diff = sasa_fluct[i] - sasa_avg[i] * sasa_avg[i] / time;

    if (diff > 0.0)
      sasa_fluctuations[i] = sqrt(diff / time);
    else
      sasa_fluctuations[i] = 0.0;
    DEBUG(10, "Average SASA of sasa atom " << i << " = " << sasa[i] << ",\tfluct = "
            << sasa_fluctuations[i]);

  }

  sasatot = sasatot_avg / time;
  const double diff = sasatot_fluct - sasatot_avg * sasatot_avg / time;
  if (diff > 0.0)
    sasatot_fluct = sqrt(diff / time);
  else
    sasatot_fluct = 0.0;
}

void configuration::Average::Block_Average
::sasavol_average(std::vector<double> & sasa_buriedvol,
                  std::vector<double> & sasa_buriedvol_fluctuations,
                  double & sasa_buriedvol_tot, double & sasa_buriedvol_totfluct)const
{
  sasa_buriedvol.resize(sasa_buriedvol_avg.size());
  sasa_buriedvol_fluctuations.resize(sasa_buriedvol_avg.size());
  
  for (unsigned int i = 0; i < sasa_buriedvol_avg.size(); ++i){
    sasa_buriedvol[i] = sasa_buriedvol_avg[i] / time;
    const double diff = sasa_buriedvol_fluct[i] - sasa_buriedvol_avg[i] * sasa_buriedvol_avg[i] / time;

    if (diff > 0.0)
      sasa_buriedvol_fluctuations[i] = sqrt(diff / time);
    else
      sasa_buriedvol_fluctuations[i] = 0.0;
    DEBUG(10, "Average volume of sasa atom " << i << " = " << sasa_buriedvol[i] <<
            "\t, fluct = " << sasa_buriedvol_fluctuations[i]);

  }

  sasa_buriedvol_tot = sasa_buriedvol_tot_avg / time;
  const double diff = sasa_buriedvol_totfluct - sasa_buriedvol_tot_avg * sasa_buriedvol_tot_avg / time;
  if (diff > 0.0)
    sasa_buriedvol_totfluct = sqrt(diff / time);
  else
    sasa_buriedvol_totfluct = 0.0;
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
  ENERGY_AVG(bond_total);
  ENERGY_AVG(angle_total);
  ENERGY_AVG(improper_total);
  ENERGY_AVG(dihedral_total);
  ENERGY_AVG(crossdihedral_total);
  ENERGY_AVG(bonded_total);
  ENERGY_AVG(nonbonded_total);
  ENERGY_AVG(lj_total);
  ENERGY_AVG(crf_total);
  ENERGY_AVG(ls_total);
  ENERGY_AVG(ls_pair_total);
  ENERGY_AVG(ls_realspace_total);
  ENERGY_AVG(ls_kspace_total);
  ENERGY_AVG(ls_self_total);
  ENERGY_AVG(ls_surface_total);
  ENERGY_AVG(ls_a_term_total);
  ENERGY_AVG(special_total);
  ENERGY_AVG(sasa_total);
  ENERGY_AVG(sasa_volume_total);
  ENERGY_AVG(qm_total);
  ENERGY_AVG(posrest_total);
  ENERGY_AVG(distanceres_total);
  ENERGY_AVG(angrest_total);
  ENERGY_AVG(dihrest_total);
  ENERGY_AVG(disfieldres_total);
  ENERGY_AVG(jvalue_total);
  ENERGY_AVG(xray_total);
  ENERGY_AVG(leus_total);
  ENERGY_AVG(oparam_total);
  ENERGY_AVG(rdc_total);   
  ENERGY_AVG(symrest_total);
  ENERGY_AVG(constraints_total);
  ENERGY_AVG(self_total);
  ENERGY_AVG(eds_vr);
  ENERGY_AVG(entropy_term);
  ENERGY_AVG(nn_valid);

  // ANITA
  for(size_t i=0; i < e.A_lj_total.size(); ++i){
    ENERGY_AVG(A_lj_total[i]); 
    ENERGY_AVG(B_lj_total[i]); 
    ENERGY_AVG(A_crf_total[i]); 
    ENERGY_AVG(B_crf_total[i]); 
    for(size_t j=0; j < e.A_lj_energy[i].size(); ++j){
      for(size_t k=0; k < e.A_lj_energy[i][j].size(); ++k){
        ENERGY_AVG(A_lj_energy[i][j][k]); 
        ENERGY_AVG(B_lj_energy[i][j][k]); 
        ENERGY_AVG(A_crf_energy[i][j][k]); 
        ENERGY_AVG(B_crf_energy[i][j][k]); 
      }
    }
  } // ANITA


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
    ENERGY_AVG(crossdihedral_energy[i]);

    // and the nonbonded groups
    assert(energy_avg.lj_energy.size() > i);
    assert(energy_fluct.lj_energy.size() > i);

    for(size_t j=0; j < e.lj_energy.size(); ++j){

      assert(energy_avg.lj_energy[i].size() > j);
      assert(energy_fluct.lj_energy[i].size() > j);

      ENERGY_AVG(lj_energy[i][j]);
      ENERGY_AVG(crf_energy[i][j]);
      ENERGY_AVG(ls_real_energy[i][j]);
      ENERGY_AVG(ls_k_energy[i][j]);
    }
  
    // special energies
    ENERGY_AVG(sasa_energy[i]);
    ENERGY_AVG(sasa_volume_energy[i]);
    ENERGY_AVG(posrest_energy[i]);
    ENERGY_AVG(distanceres_energy[i]);
    ENERGY_AVG(angrest_energy[i]);
    ENERGY_AVG(dihrest_energy[i]);
    ENERGY_AVG(disfieldres_energy[i]);
    ENERGY_AVG(constraints_energy[i]);
    ENERGY_AVG(jvalue_energy[i]);
    ENERGY_AVG(rdc_energy[i]);
    ENERGY_AVG(self_energy[i]);
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

  // pressure
  pressure_avg = old.pressure_avg + conf.old().pressure_tensor * dt;
  pressure_fluct = old.pressure_fluct + math::square(conf.old().pressure_tensor) * dt;
  // virial
  virial_avg = old.virial_avg + conf.old().virial_tensor * dt;
  virial_fluct = old.virial_fluct + math::square(conf.old().virial_tensor) * dt;
  // ekin_tensor
  ekin_tensor_avg = old.ekin_tensor_avg + conf.old().kinetic_energy_tensor * dt;
  ekin_tensor_fluct = old.ekin_tensor_fluct + math::square(conf.old().kinetic_energy_tensor) * dt;

  // box
  box_avg = old.box_avg + conf.old().box * dt;
  box_fluct = old.box_fluct + math::square(conf.old().box) * dt;

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
::update_sasa(topology::Topology const & topo,
              configuration::Configuration const & conf,
              simulation::Simulation const & sim,
	      Block_Average const & old)
{
  const double dt = sim.time_step_size();

  // average sasa and fluctuation for each atom i
  const unsigned int num_sasa_atoms = topo.sasa_parameter().size();
  for (unsigned int i = 0; i < num_sasa_atoms; ++i){
    sasa_avg[i] = old.sasa_avg[i] + dt * conf.old().sasa_area[i];
    sasa_fluct[i] = old.sasa_fluct[i] + conf.old().sasa_area[i] * conf.old().sasa_area[i] * dt;
  }
  // average total sasa and fluctuation
  sasatot_avg = old.sasatot_avg + conf.old().sasa_tot * dt;
  sasatot_fluct = old.sasatot_fluct + conf.old().sasa_tot * conf.old().sasa_tot * dt;
  
  DEBUG(10, "Updated average total SASA: " << sasatot_avg <<
          "\tand old average total SASA: " << old.sasatot_avg << "\tat time: " << time);
}

void configuration::Average::Block_Average
::update_sasavol(topology::Topology const & topo,
                 configuration::Configuration const & conf,
                 simulation::Simulation const & sim,
	         Block_Average const & old)
{
  const double dt = sim.time_step_size();

  // average volume and fluctuation for atom i 
  const unsigned int num_sasa_atoms = topo.sasa_parameter().size();
  for (unsigned int i = 0; i < num_sasa_atoms; ++i){
    sasa_buriedvol_avg[i] = old.sasa_buriedvol_avg[i] + dt * conf.old().sasa_buriedvol[i];
    sasa_buriedvol_fluct[i] = old.sasa_buriedvol_fluct[i] + conf.old().sasa_buriedvol[i] * conf.old().sasa_buriedvol[i] * dt;
  }
  // average total volume and fluctuation
  sasa_buriedvol_tot_avg = old.sasa_buriedvol_tot_avg + conf.old().sasa_buriedvol_tot * dt;
  sasa_buriedvol_tot_fluct = old.sasa_buriedvol_tot_fluct + conf.old().sasa_buriedvol_tot * conf.old().sasa_buriedvol_tot * dt;

  DEBUG(10, "Updated average total volume: " << sasa_buriedvol_tot_avg <<
          "\tand old average total volume: " << old.sasa_buriedvol_tot_avg << "\tat time: " << time);
}

void configuration::Average::Block_Average
::energy_average_helper(configuration::Energy const & avg,
			configuration::Energy const & fluct,
			configuration::Energy &energy, 
			configuration::Energy &fluctuation,
			double const dlamt)const
{
  DEBUG(10, "new average total: " << energy_avg.total << "\ttime: " << time);

  double diff = 0.0;
  double cumu = time;
  if (dlamt) cumu = 1;
  
  Energy &e = energy;
  Energy &f = fluctuation;

  DEBUG(5, "ANITA energy_average_helper, resizing to " <<
            energy_avg.A_lj_total.size());
  e.resize(unsigned(energy_avg.bond_energy.size()),
		unsigned(energy_avg.kinetic_energy.size()),
                unsigned(energy_avg.A_lj_total.size())); //ANITA
  f.resize(unsigned(energy_avg.bond_energy.size()),
	  unsigned(energy_avg.kinetic_energy.size()),
          unsigned(energy_avg.A_lj_total.size())); //ANITA

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
  ENERGY_RES(crossdihedral_total);

  ENERGY_RES(nonbonded_total);
  ENERGY_RES(lj_total);
  ENERGY_RES(crf_total);
  ENERGY_RES(ls_total);
  ENERGY_RES(ls_pair_total);
  ENERGY_RES(ls_realspace_total);
  ENERGY_RES(ls_kspace_total);
  ENERGY_RES(ls_self_total);
  ENERGY_RES(ls_surface_total);
  ENERGY_RES(ls_a_term_total);

  ENERGY_RES(special_total);
  ENERGY_RES(sasa_total);
  ENERGY_RES(sasa_volume_total);  
  ENERGY_RES(qm_total);
  ENERGY_RES(posrest_total);
  ENERGY_RES(distanceres_total);
  ENERGY_RES(angrest_total);
  ENERGY_RES(dihrest_total);
  ENERGY_RES(disfieldres_total);
  ENERGY_RES(jvalue_total);
  ENERGY_RES(xray_total);
  ENERGY_RES(leus_total);
  ENERGY_RES(oparam_total);
  ENERGY_RES(rdc_total); 
  ENERGY_RES(symrest_total);
  ENERGY_RES(constraints_total);
  ENERGY_RES(entropy_term);
  ENERGY_RES(self_total);
  ENERGY_RES(nn_valid);

  // ANITA
  for(size_t i=0; i < e.A_lj_total.size(); ++i){
    ENERGY_RES(A_lj_total[i]);
    ENERGY_RES(B_lj_total[i]);
    ENERGY_RES(A_crf_total[i]);
    ENERGY_RES(B_crf_total[i]);
    for(size_t j=0; j < e.A_lj_energy[i].size(); ++j){
      for(size_t k=0; k < e.A_lj_energy[i][j].size(); ++k){
        ENERGY_RES(A_lj_energy[i][j][k]);
        ENERGY_RES(B_lj_energy[i][j][k]);
        ENERGY_RES(A_crf_energy[i][j][k]);
        ENERGY_RES(B_crf_energy[i][j][k]);
      }
    }
  } // ANITA

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
    ENERGY_RES(crossdihedral_energy[i]);

    // and the nonbonded groups
    assert(energy_avg.lj_energy.size() > i);
    assert(energy_fluct.lj_energy.size() > i);

    for(size_t j=0; j < e.lj_energy.size(); ++j){

      assert(energy_avg.lj_energy[i].size() > j);
      assert(energy_fluct.lj_energy[i].size() > j);

      ENERGY_RES(lj_energy[i][j]);
      ENERGY_RES(crf_energy[i][j]);
      ENERGY_RES(ls_real_energy[i][j]);
      ENERGY_RES(ls_k_energy[i][j]);
    }
  
    // special energies
    ENERGY_RES(sasa_energy[i]);
    ENERGY_RES(sasa_volume_energy[i]);
    ENERGY_RES(posrest_energy[i]);
    ENERGY_RES(distanceres_energy[i]);
    ENERGY_RES(angrest_energy[i]);
    ENERGY_RES(dihrest_energy[i]);
    ENERGY_RES(disfieldres_energy[i]);
    ENERGY_RES(constraints_energy[i]);
    ENERGY_RES(jvalue_energy[i]);
    ENERGY_RES(rdc_energy[i]);
    
    ENERGY_RES(self_energy[i]);
  }

#undef ENERGY_RES
}
