/**
 * @file energy_average.tcc
 * implements member methods.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE energy

#include <util/stdheader.h>
#include <configuration/configuration_global.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>

configuration::Energy_Average::Energy_Average()
  : m_average(),
    m_square_average(),
    m_time(0)
{
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      m_pressure_average(a,b) = 0.0;
      m_square_pressure_average(a,b) = 0.0;
    }
  }
  
}

void
configuration::Energy_Average::zero()
{
  DEBUG(10, "energy average: zero");
  m_average.zero();
  m_square_average.zero();
  m_time = 0.0;

  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      m_pressure_average(a,b) = 0.0;
      m_square_pressure_average(a,b) = 0.0;
    }
  }

}

void
configuration::Energy_Average::resize(size_t const energy_groups, size_t const multi_baths)
{
  DEBUG(10, "energy average: resize");
  m_average.resize(energy_groups, multi_baths);
  m_square_average.resize(energy_groups, multi_baths);
  m_time = 0.0;
}

void
configuration::Energy_Average::update(configuration::Energy const &e, 
				      configuration::Energy_Average const &av,
				      double const dt)
{
  // totals
  DEBUG(10, "energy average: update - kin size: " << e.kinetic_energy.size());
  
  DEBUG(10, "old total: " << av.m_average.total
	<< "\ttotal: " << e.total << "\tdt: " << dt);
  
  m_average.total = av.m_average.total + dt * e.total;
  
  m_square_average.total = av.m_square_average.total + dt * e.total * e.total;
  
  m_average.kinetic_total = av.m_average.kinetic_total + dt * e.kinetic_total;
  m_square_average.kinetic_total = av.m_square_average.kinetic_total + 
    dt * e.kinetic_total * e.kinetic_total;
  
  m_average.potential_total = av.m_average.potential_total +
    dt * e.potential_total;
  m_square_average.potential_total = av.m_square_average.potential_total + 
    dt * e.potential_total * e.potential_total;

  // physical interaction totals
  m_average.bond_total = av.m_average.bond_total + dt * e.bond_total;
  m_square_average.bond_total = av.m_square_average.bond_total + 
    dt * e.bond_total * e.bond_total;
  
  m_average.angle_total = av.m_average.angle_total + dt * e.angle_total;
  m_square_average.angle_total = av.m_square_average.angle_total +
    dt * e.angle_total * e.angle_total;
  
  m_average.improper_total = av.m_average.improper_total +
    dt * e.improper_total;
  m_square_average.improper_total = av.m_square_average.improper_total + 
    dt * e.improper_total * e.improper_total;
  
  m_average.dihedral_total = av.m_average.dihedral_total + dt * e.dihedral_total;
  m_square_average.dihedral_total = av.m_square_average.dihedral_total +
    dt * e.dihedral_total * e.dihedral_total;
  
  m_average.bonded_total = av.m_average.bonded_total + dt * e.bonded_total;
  m_square_average.bonded_total = av.m_square_average.bonded_total +
    dt * e.bonded_total * e.bonded_total;
  
  m_average.nonbonded_total = av.m_average.nonbonded_total +
    dt * e.nonbonded_total;
  m_square_average.nonbonded_total = av.m_square_average.nonbonded_total +
    dt * e.nonbonded_total * e.nonbonded_total;
  
  m_average.lj_total = av.m_average.lj_total + dt * e.lj_total;
  m_square_average.lj_total = av.m_square_average.lj_total +
    dt * e.lj_total * e.lj_total;
  
  m_average.crf_total = av.m_average.crf_total + dt * e.crf_total;
  m_square_average.crf_total = av.m_square_average.crf_total +
    dt * e.crf_total * e.crf_total;
  
  m_average.special_total = av.m_average.special_total + dt * e.special_total;
  m_square_average.special_total = av.m_square_average.special_total +
    dt * e.special_total * e.special_total;

  m_average.posrest_total = av.m_average.posrest_total + dt * e.posrest_total;
  m_square_average.posrest_total = av.m_square_average.posrest_total +
    dt * e.posrest_total * e.posrest_total;
  
  // kinetic energies of the baths
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){

    assert(m_average.kinetic_energy.size() > i);
    assert(m_square_average.kinetic_energy.size() > i);
    assert(e.kinetic_energy.size() > i);
    
    m_average.kinetic_energy[i] = av.m_average.kinetic_energy[i] +
      dt * e.kinetic_energy[i];
    m_square_average.kinetic_energy[i] = av.m_square_average.kinetic_energy[i] + 
      dt * e.kinetic_energy[i] * e.kinetic_energy[i];

    assert(m_average.com_kinetic_energy.size() > i);
    assert(m_square_average.com_kinetic_energy.size() > i);
    assert(e.com_kinetic_energy.size() > i);
    
    m_average.com_kinetic_energy[i] = av.m_average.com_kinetic_energy[i] +
      dt * e.com_kinetic_energy[i];
    m_square_average.com_kinetic_energy[i] = av.m_square_average.com_kinetic_energy[i] +
      dt * e.com_kinetic_energy[i] * e.com_kinetic_energy[i];

    assert(m_average.ir_kinetic_energy.size() > i);
    assert(m_square_average.ir_kinetic_energy.size() > i);
    assert(e.ir_kinetic_energy.size() > i);

    m_average.ir_kinetic_energy[i] = av.m_average.ir_kinetic_energy[i] +
      dt * e.ir_kinetic_energy[i];
    m_square_average.ir_kinetic_energy[i] = av.m_square_average.ir_kinetic_energy[i] +
      dt * e.ir_kinetic_energy[i] * e.ir_kinetic_energy[i];
    
  }
  
  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){
    m_average.bond_energy[i] = av.m_average.bond_energy[i] + dt * e.bond_energy[i];
    m_square_average.bond_energy[i] = av.m_square_average.bond_energy[i] +
      dt * e.bond_energy[i] * e.bond_energy[i];
  
    m_average.angle_energy[i] = av.m_average.angle_energy[i] + dt * e.angle_energy[i];
    m_square_average.angle_energy[i] = av.m_square_average.angle_energy[i] +
      dt * e.angle_energy[i] * e.angle_energy[i];
  
    m_average.improper_energy[i] = av.m_average.improper_energy[i] +
      dt * e.improper_energy[i];
    m_square_average.improper_energy[i] = av.m_square_average.improper_energy[i] +
      dt * e.improper_energy[i] * e.improper_energy[i];
  
    m_average.dihedral_energy[i] = av.m_average.dihedral_energy[i] +
      dt * e.dihedral_energy[i];
    m_square_average.dihedral_energy[i] = av.m_square_average.dihedral_energy[i] +
      dt * e.dihedral_energy[i] * e.dihedral_energy[i];

    // and the nonbonded groups
    for(size_t j=0; j < e.lj_energy.size(); ++j){
      m_average.lj_energy[i][j] = av.m_average.lj_energy[i][j] + 
	dt * e.lj_energy[i][j];
      m_square_average.lj_energy[i][j] = av.m_square_average.lj_energy[i][j] +
	dt * e.lj_energy[i][j] * e.lj_energy[i][j];
  
      m_average.crf_energy[i][j] = av.m_average.crf_energy[i][j] +
	dt * e.crf_energy[i][j];
      m_square_average.crf_energy[i][j] = av.m_square_average.crf_energy[i][j] +
	dt * e.crf_energy[i][j] * e.crf_energy[i][j];
      
    }
  
    // special energies
    m_average.posrest_energy[i] = av.m_average.posrest_energy[i] + dt * e.posrest_energy[i];
    m_square_average.posrest_energy[i] = av.m_square_average.posrest_energy[i] +
      dt * e.posrest_energy[i] * e.posrest_energy[i];
   

  }

  m_time = av.m_time + dt;
  
}

void
configuration::Energy_Average::update(math::Matrix const & pressure,
				      configuration::Energy_Average const &av,
				      double const dt)
{
  // the pressure...
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      m_pressure_average(a,b) = av.m_pressure_average(a,b) + dt * pressure(a,b);
      m_square_pressure_average(a,b) = av.m_square_pressure_average(a,b) +
	dt * pressure(a,b) * pressure(a,b);
    }
  }
}

void
configuration::Energy_Average::average(configuration::Energy &energy, 
				       configuration::Energy &fluctuation,
				       math::Matrix &pressure,
				       math::Matrix &pressure_fluctuations)
{
  double diff;
  
  Energy &e = energy;
  Energy &f = fluctuation;

  e.resize(m_average.bond_energy.size(), m_average.kinetic_energy.size());
  f.resize(m_average.bond_energy.size(), m_average.kinetic_energy.size());

  DEBUG(10, "average total: " << m_average.total << "\ttime: " << m_time);

  e.total = m_average.total / m_time;
  diff = m_square_average.total - m_average.total * m_average.total / m_time;
  if (diff > 0.0)
    f.total = sqrt(diff / m_time);
  else f.total = 0.0;
  
  e.kinetic_total = m_average.kinetic_total / m_time;
  diff = m_square_average.kinetic_total -
    m_average.kinetic_total * m_average.kinetic_total / m_time;
  if (diff > 0.0)
    f.kinetic_total = sqrt(diff / m_time);
  else
    f.kinetic_total = 0.0;
  
  e.potential_total = m_average.potential_total / m_time;
  diff = m_square_average.potential_total - 
    m_average.potential_total * m_average.potential_total / m_time;
  if (diff > 0.0)
    f.potential_total = sqrt(diff / m_time);
  else
    f.potential_total = 0.0;
  
  // physical interaction totals
  e.bond_total = m_average.bond_total / m_time;
  diff = m_square_average.bond_total -
    m_average.bond_total * m_average.bond_total / m_time;
  if (diff > 0.0)
    f.bond_total = sqrt(diff / m_time);
  else
    f.bond_total = 0.0;
  
  e.angle_total = m_average.angle_total / m_time;
  diff = m_square_average.angle_total -
    m_average.angle_total * m_average.angle_total / m_time;
  if (diff > 0.0)
    f.angle_total = sqrt(diff / m_time);
  else
    f.angle_total = 0.0;
  
  e.improper_total = m_average.improper_total / m_time;
  diff = m_square_average.improper_total -
    m_average.improper_total * m_average.improper_total / m_time;
  if (diff > 0.0)
    f.improper_total = sqrt(diff / m_time);
  else
    f.improper_total = 0.0;
  
  e.dihedral_total = m_average.dihedral_total / m_time;
  diff = m_square_average.dihedral_total -
    m_average.dihedral_total * m_average.dihedral_total / m_time;
  if (diff > 0.0)
    f.dihedral_total = sqrt(diff / m_time);
  else
    f.dihedral_total = 0;

  e.bonded_total = m_average.bonded_total / m_time;
  diff = m_square_average.bonded_total -
    m_average.bonded_total * m_average.bonded_total / m_time;
  if (diff > 0.0)
    f.bonded_total = sqrt(diff / m_time);
  else
    f.bonded_total = 0.0;

  e.nonbonded_total = m_average.nonbonded_total / m_time;
  diff = m_square_average.nonbonded_total -
    m_average.nonbonded_total * m_average.nonbonded_total / m_time;
  if (diff > 0.0)
    f.nonbonded_total = sqrt(diff / m_time);
  else
    f.nonbonded_total = 0.0;
  
  e.lj_total = m_average.lj_total / m_time;
  diff = m_square_average.lj_total -
    m_average.lj_total * m_average.lj_total / m_time;
  if (diff > 0.0)
    f.lj_total = sqrt(diff / m_time);
  else
    f.lj_total = 0.0;

  e.crf_total = m_average.crf_total / m_time;
  diff = m_square_average.crf_total -
    m_average.crf_total * m_average.crf_total / m_time;
  if (diff > 0.0)
    f.crf_total = sqrt(diff / m_time);
  else
    f.crf_total = 0.0;

  // special energies
  e.special_total = m_average.special_total / m_time;
  diff = m_square_average.special_total -
    m_average.special_total * m_average.special_total / m_time;
  if (diff > 0.0)
    f.special_total = sqrt(diff / m_time);
  else
    f.special_total = 0.0;

  e.posrest_total = m_average.posrest_total / m_time;
  diff = m_square_average.posrest_total -
    m_average.posrest_total * m_average.posrest_total / m_time;
  if (diff > 0.0)
    f.posrest_total = sqrt(diff / m_time);
  else
    f.posrest_total = 0.0;
  
  // kinetic energies of the baths
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){
    e.kinetic_energy[i] = m_average.kinetic_energy[i] / m_time;
    diff = m_square_average.kinetic_energy[i] -
      m_average.kinetic_energy[i] * m_average.kinetic_energy[i] / m_time;
    if (diff > 0.0)
      f.kinetic_energy[i] = sqrt(diff / m_time);
    else
      f.com_kinetic_energy[i] = 0.0;

    e.com_kinetic_energy[i] = m_average.com_kinetic_energy[i] / m_time;
    diff = m_square_average.com_kinetic_energy[i] -
      m_average.com_kinetic_energy[i] * m_average.com_kinetic_energy[i] 
      / m_time;
    if (diff > 0.0)
      f.com_kinetic_energy[i] = sqrt(diff / m_time);
    else
      f.com_kinetic_energy[i] = 0.0;

    e.ir_kinetic_energy[i] = m_average.ir_kinetic_energy[i] / m_time;
    diff = m_square_average.ir_kinetic_energy[i] -
      m_average.ir_kinetic_energy[i] * m_average.ir_kinetic_energy[i] 
      / m_time;
    if (diff > 0.0)
      f.ir_kinetic_energy[i] = sqrt(diff / m_time);
    else
      f.ir_kinetic_energy[i] = 0.0;
    
  }
  
  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){
    e.bond_energy[i] = m_average.bond_energy[i] / m_time;
    diff = m_square_average.bond_energy[i] -
      m_average.bond_energy[i] * m_average.bond_energy[i] / m_time;
    if (diff > 0.0)
      f.bond_energy[i] = sqrt(diff / m_time);
    else
      f.bond_energy[i] = 0.0;

    e.angle_energy[i] = m_average.angle_energy[i] / m_time;
    diff = m_square_average.angle_energy[i] -
      m_average.angle_energy[i] * m_average.angle_energy[i] / m_time;
    if (diff > 0.0)
      f.angle_energy[i] = sqrt(diff / m_time);
    else
      f.angle_energy[i] = 0.0;

    e.improper_energy[i] = m_average.improper_energy[i] / m_time;
    diff = m_square_average.improper_energy[i]
      - m_average.improper_energy[i] * m_average.improper_energy[i] / m_time;
    if (diff > 0.0)
      f.improper_energy[i] = sqrt(diff / m_time);
    else
      f.improper_energy[i] = 0.0;
    
    e.dihedral_energy[i] = m_average.dihedral_energy[i] / m_time;
    diff = m_square_average.dihedral_energy[i] -
      m_average.dihedral_energy[i] * m_average.dihedral_energy[i] / m_time;
    if (diff > 0.0)
      f.dihedral_energy[i] = sqrt(diff / m_time);
    else
      f.dihedral_energy[i] = 0.0;
    
    // and the nonbonded groups
    for(size_t j=0; j < e.lj_energy.size(); ++j){
      e.lj_energy[i][j] = m_average.lj_energy[i][j] / m_time;
      diff = m_square_average.lj_energy[i][j] -
	m_average.lj_energy[i][j] * m_average.lj_energy[i][j] / m_time;
      if (diff > 0.0)
	f.lj_energy[i][j] = sqrt(diff / m_time);
      else
	f.lj_energy[i][j] = 0.0;

      e.crf_energy[i][j] = m_average.crf_energy[i][j] / m_time;
      diff = m_square_average.crf_energy[i][j] -
	m_average.crf_energy[i][j] * m_average.crf_energy[i][j] / m_time;
      if (diff > 0.0)
	f.crf_energy[i][j] = sqrt(diff / m_time);
      else
	f.crf_energy[i][j] = 0.0;
    }

    e.posrest_energy[i] = m_average.posrest_energy[i] / m_time;
    diff = m_square_average.posrest_energy[i]
      - m_average.posrest_energy[i] * m_average.posrest_energy[i] / m_time;
    if (diff > 0.0)
      f.posrest_energy[i] = sqrt(diff / m_time);
    else
      f.posrest_energy[i] = 0.0;

  
  }

  // the pressure
  for(int a=0; a<3; ++a){
    for(int b=0; b<3; ++b){
      pressure(a,b) = m_pressure_average(a,b) / m_time;
      diff = m_square_pressure_average(a,b) -
	m_pressure_average(a,b) * m_pressure_average(a,b) / m_time;
      if (diff > 0.0)
	pressure_fluctuations(a,b) = sqrt(diff / m_time);
      else
	pressure_fluctuations(a,b) = 0.0;
    }
  }

}
