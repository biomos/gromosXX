/**
 * @file energy.tcc
 * implements the energy methods.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE energy

#include <util/stdheader.h>
#include <configuration/configuration_global.h>
#include <configuration/energy.h>

void configuration::Energy::zero(bool potential, bool kinetic)
{
  DEBUG(10, "energy: zero");
  
  if (potential){
    DEBUG(15, "zero potential energies");

    total = 0.0;
    potential_total = 0.0;
    bond_total = 0.0;
    angle_total = 0.0;
    improper_total = 0.0;
    dihedral_total = 0.0;
    bonded_total = 0.0;
    nonbonded_total = 0.0;
    lj_total = 0.0;
    crf_total = 0.0;
    special_total = 0.0;
    posrest_total = 0.0;
    
    bond_energy.assign(bond_energy.size(), 0.0);
    angle_energy.assign(angle_energy.size(), 0.0);
    improper_energy.assign(improper_energy.size(), 0.0);
    dihedral_energy.assign(dihedral_energy.size(), 0.0);
    posrest_energy.assign(posrest_energy.size(), 0.0);

    DEBUG(15, "energy groups: " << lj_energy.size() << " - " << crf_energy.size());

    lj_energy.assign(lj_energy.size(), 
		     std::vector<double>(lj_energy.size(), 0.0));
    
    crf_energy.assign(crf_energy.size(), 
		      std::vector<double>(crf_energy.size(), 0.0));
    
  }

  if (kinetic){

    DEBUG(15, "zero kinetic energies");
    
    kinetic_total = 0.0;

    kinetic_energy.assign(kinetic_energy.size(), 0.0);
    com_kinetic_energy.assign(com_kinetic_energy.size(), 0.0);
    ir_kinetic_energy.assign(ir_kinetic_energy.size(), 0.0);

  }
  DEBUG(15, "energies zero");

}


void configuration::Energy::resize(size_t const energy_groups, size_t const multi_baths)
{
  DEBUG(10, "energy resize");
  
  if (energy_groups){
    bond_energy.resize(energy_groups);
    angle_energy.resize(energy_groups);
    improper_energy.resize(energy_groups);
    dihedral_energy.resize(energy_groups);
  
    lj_energy.resize(energy_groups);
    crf_energy.resize(energy_groups);

    posrest_energy.resize(energy_groups);
    
    for(size_t i=0; i<energy_groups; ++i){
      lj_energy[i].resize(energy_groups);
      crf_energy[i].resize(energy_groups);  
    }
  }

  if (multi_baths){
    kinetic_energy.resize(multi_baths);
    com_kinetic_energy.resize(multi_baths);
    ir_kinetic_energy.resize(multi_baths);
  }

  zero();  
}

void configuration::Energy::calculate_totals()
{
  DEBUG(10, "energy: calculate totals");
  
  int num_groups = bond_energy.size();

  kinetic_total = 0.0;

  bond_total = 0.0;
  angle_total = 0.0;
  improper_total = 0.0;
  dihedral_total = 0.0;
  lj_total = 0.0;
  crf_total = 0.0;
    
  for(size_t i=0; i<kinetic_energy.size(); ++i){
    kinetic_total += kinetic_energy[i];
  }


  for(int i=0; i<num_groups; i++){
    for(int j=0; j<num_groups; j++){
      lj_total   += lj_energy[i][j];
      crf_total  += crf_energy[i][j];
    }

    bond_total     += bond_energy[i];
    angle_total    += angle_energy[i];
    improper_total += improper_energy[i];
    dihedral_total += dihedral_energy[i];
    posrest_total += posrest_energy[i];
  }

  nonbonded_total = lj_total + crf_total;
  bonded_total    = bond_total + angle_total + 
    dihedral_total + improper_total;
  potential_total = nonbonded_total + bonded_total;
  
  special_total = posrest_total;

  total = potential_total + kinetic_total + special_total;

}
