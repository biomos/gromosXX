/**
 * @file energy.cc
 * implements the energy methods.
 */

#include "../stdheader.h"

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE energy

#include "../configuration/configuration_global.h"
#include "../configuration/energy.h"

#include "../util/error.h"

#include "energy.h"
#include <vector>

configuration::Energy::Energy() :
total(0.0),
kinetic_total(0.0),
potential_total(0.0),
bond_total(0.0),
angle_total(0.0),
improper_total(0.0),
dihedral_total(0.0),
crossdihedral_total(0.0),
bonded_total(0.0),
nonbonded_total(0.0),
lj_total(0.0),
crf_total(0.0),
ls_total(0.0),
ls_pair_total(0.0),
ls_realspace_total(0.0),
ls_kspace_total(0.0),
ls_self_total(0.0),
ls_surface_total(0.0),
ls_a_term_total(0.0),
special_total(0.0),
posrest_total(0.0),
distanceres_total(0.0),
dihrest_total(0.0),
jvalue_total(0.0),
xray_total(0.0),
leus_total(0.0),
constraints_total(0.0),
external_total(0.0),
self_total(0.0),
sasa_total(0.0),
sasa_volume_total(0.0),
qm_total(0.0),
eds_vr(0.0),
entropy_term(0.0),
m_ewarn(1E99) {
}

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
    crossdihedral_total = 0.0;
    bonded_total = 0.0;
    nonbonded_total = 0.0;
    lj_total = 0.0;
    crf_total = 0.0;
    ls_total = 0.0;
    ls_pair_total = 0.0;
    ls_realspace_total = 0.0;
    ls_kspace_total = 0.0;
    ls_self_total = 0.0;
    ls_surface_total = 0.0;
    ls_a_term_total = 0.0;
    special_total = 0.0;
    posrest_total = 0.0;
    distanceres_total = 0.0;
    dihrest_total = 0.0;
    jvalue_total = 0.0;
    xray_total = 0.0;
    leus_total = 0.0;
    constraints_total = 0.0;
    entropy_term = 0.0;
    external_total = 0.0;
    self_total = 0.0;
    sasa_total = 0.0;
    sasa_volume_total = 0.0;
    qm_total = 0.0;
    eds_vr = 0.0;
    eds_vi.assign(eds_vi.size(), 0.0);
    eds_vi_special.assign(eds_vi_special.size(), 0.0);
    
    bond_energy.assign(bond_energy.size(), 0.0);
    angle_energy.assign(angle_energy.size(), 0.0);
    improper_energy.assign(improper_energy.size(), 0.0);
    dihedral_energy.assign(dihedral_energy.size(), 0.0);
    crossdihedral_energy.assign(crossdihedral_energy.size(), 0.0);
    posrest_energy.assign(posrest_energy.size(), 0.0);
    distanceres_energy.assign(distanceres_energy.size(), 0.0);
    dihrest_energy.assign(dihrest_energy.size(), 0.0);
    jvalue_energy.assign(jvalue_energy.size(), 0.0);
    constraints_energy.assign(constraints_energy.size(), 0.0);
    self_energy.assign(self_energy.size(), 0.0);
    sasa_energy.assign(sasa_energy.size(), 0.0);
    sasa_volume_energy.assign(sasa_volume_energy.size(), 0.0);

    DEBUG(15, "energy groups: " << unsigned(lj_energy.size()) 
			<< " - " << unsigned(crf_energy.size()));

    lj_energy.assign(lj_energy.size(), 
		     std::vector<double>(lj_energy.size(), 0.0));
    
    crf_energy.assign(crf_energy.size(), 
		      std::vector<double>(crf_energy.size(), 0.0));
    
    ls_real_energy.assign(ls_real_energy.size(),
            std::vector<double>(ls_real_energy.size(), 0.0));
    
    ls_k_energy.assign(ls_k_energy.size(),
            std::vector<double>(ls_k_energy.size(), 0.0));
    
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


void configuration::Energy::resize(unsigned int energy_groups, unsigned int multi_baths)
{
  DEBUG(10, "energy resize");
  
  if (energy_groups){
    bond_energy.resize(energy_groups);
    angle_energy.resize(energy_groups);
    improper_energy.resize(energy_groups);
    dihedral_energy.resize(energy_groups);
    crossdihedral_energy.resize(energy_groups);
  
    lj_energy.resize(energy_groups);
    crf_energy.resize(energy_groups);
    ls_real_energy.resize(energy_groups);
    ls_k_energy.resize(energy_groups);

    posrest_energy.resize(energy_groups);
    distanceres_energy.resize(energy_groups);
    dihrest_energy.resize(energy_groups);
    jvalue_energy.resize(energy_groups);
    constraints_energy.resize(energy_groups);
    
    self_energy.resize(energy_groups);

    sasa_energy.resize(energy_groups);
    sasa_volume_energy.resize(energy_groups);
    
    for(unsigned int i=0; i<energy_groups; ++i){
      lj_energy[i].resize(energy_groups);
      crf_energy[i].resize(energy_groups);
      ls_real_energy[i].resize(energy_groups);
      ls_k_energy[i].resize(energy_groups);
    }
  }

  if (multi_baths){
    kinetic_energy.resize(multi_baths);
    com_kinetic_energy.resize(multi_baths);
    ir_kinetic_energy.resize(multi_baths);
  }

  zero();  
}

int configuration::Energy::calculate_totals()
{
  DEBUG(10, "energy: calculate totals");
  
  int num_groups = unsigned(bond_energy.size());

  kinetic_total = 0.0;

  bond_total = 0.0;
  angle_total = 0.0;
  improper_total = 0.0;
  dihedral_total = 0.0;
  crossdihedral_total = 0.0;
  lj_total = 0.0;
  crf_total = 0.0;
  ls_total = 0.0;
  ls_realspace_total = 0.0;
  //ls_kspace_total = 0.0;
  posrest_total = 0.0; 
  distanceres_total =0.0; 
  dihrest_total = 0.0;
  jvalue_total = 0.0; 
  constraints_total = 0.0;
  self_total = 0.0;
  sasa_total = 0.0;
  sasa_volume_total = 0.0;
  
  for(size_t i=0; i<kinetic_energy.size(); ++i){
    if (kinetic_energy[i] > m_ewarn){
      std::cout << "EWARN: kinetic energy " << i+1 << " = " << kinetic_energy[i] << "\n";
    }
    
    kinetic_total += kinetic_energy[i];
  }


  for(int i=0; i<num_groups; i++){
    for(int j=i; j<num_groups; j++){

      if(i!=j){
        lj_energy[j][i] += lj_energy[i][j];
        crf_energy[j][i] += crf_energy[i][j];
        ls_real_energy[j][i] += ls_real_energy[i][j];
        ls_k_energy[j][i] += ls_k_energy[i][j];
        
        lj_energy[i][j] = 0.0;
        crf_energy[i][j] = 0.0;
        ls_real_energy[i][j] = 0.0;
        ls_k_energy[i][j] = 0.0;
      }

      if (lj_energy[i][j] > m_ewarn){
        std::cout << "EWARN: lj energy " << i+1 << ", " << j+1 << " = " << lj_energy[j][i] << "\n";
      }
      if (crf_energy[i][j] > m_ewarn){
        std::cout << "EWARN: crf energy " << i+1 << ", " << j+1 << " = " << crf_energy[j][i] << "\n";
      }
      if (ls_real_energy[i][j] > m_ewarn){
        std::cout << "EWARN: lj energy " << i+1 << ", " << j+1 << " = " << ls_real_energy[j][i] << "\n";
      }
      if (ls_k_energy[i][j] > m_ewarn){
        std::cout << "EWARN: crf energy " << i+1 << ", " << j+1 << " = " << ls_k_energy[j][i] << "\n";
      }
      
      lj_total   += lj_energy[j][i];
      crf_total  += crf_energy[j][i];
      ls_realspace_total   += ls_real_energy[j][i];
      //ls_kspace_total  += ls_k_energy[j][i];
    }

    if (bond_energy[i] > m_ewarn){
      std::cout << "EWARN: bond energy " << i+1 << " = " << bond_energy[i] << "\n";
    }
    bond_total         += bond_energy[i];
    if (angle_energy[i] > m_ewarn){
      std::cout << "EWARN: angle energy " << i+1 << " = " << angle_energy[i] << "\n";
    }
    angle_total        += angle_energy[i];
    if (improper_energy[i] > m_ewarn){
      std::cout << "EWARN: improper energy " << i+1 << " = " << improper_energy[i] << "\n";
    }
    improper_total     += improper_energy[i];
    if (dihedral_energy[i] > m_ewarn){
      std::cout << "EWARN: dihedral energy " << i+1 << " = " << dihedral_energy[i] << "\n";
    }
    dihedral_total     += dihedral_energy[i];
    if (crossdihedral_energy[i] > m_ewarn){
      std::cout << "EWARN: crossdihedral energy " << i+1 << " = " << crossdihedral_energy[i] << "\n";
    }
    crossdihedral_total     += crossdihedral_energy[i];
    if (posrest_energy[i] > m_ewarn){
      std::cout << "EWARN: posrest energy " << i+1 << " = " << posrest_energy[i] << "\n";
    }
    posrest_total      += posrest_energy[i];
    if (distanceres_energy[i] > m_ewarn){
      std::cout << "EWARN: distanceres energy " << i+1 << " = " << distanceres_energy[i] << "\n";
    }
    distanceres_total     += distanceres_energy[i];
    if (dihrest_energy[i] > m_ewarn){
      std::cout << "EWARN: dihrest energy " << i+1 << " = " << dihrest_energy[i] << "\n";
    }
    dihrest_total      += dihrest_energy[i];
    if (jvalue_energy[i] > m_ewarn){
      std::cout << "EWARN: jvalue energy " << i+1 << " = " << jvalue_energy[i] << "\n";
    }
    jvalue_total       += jvalue_energy[i];
    if (constraints_energy[i] > m_ewarn){
      std::cout << "EWARN: constraints energy " << i+1 << " = " << constraints_energy[i] << "\n";
    }
    constraints_total  += constraints_energy[i];
    if (self_energy[i] > m_ewarn){
       std::cout << "EWARN: self energy " << i+1 << " = " << self_energy[i] << "\n";
    }
    self_total += self_energy[i];
    if (sasa_energy[i] > m_ewarn){
       std::cout << "EWARN: sasa energy " << i+1 << " = " << sasa_energy[i] << "\n";
    }
    sasa_total += sasa_energy[i];
    if (sasa_volume_energy[i] > m_ewarn) {
      std::cout << "EWARN: sasa volume energy " << i + 1 << " = " << sasa_volume_energy[i] << "\n";
    }
    sasa_volume_total += sasa_volume_energy[i];
  }

  if (xray_total > m_ewarn) {
    std::cout << "EWARN: xray energy = " << xray_total << "\n";
  }

  if (leus_total > m_ewarn) {
    std::cout << "EWARN: local elevation energy = " << leus_total << "\n";
  }

  ls_pair_total = ls_realspace_total + ls_kspace_total + ls_a_term_total;
  //        E(pair) + DeltaG(self) + DeltaG(surf)
  ls_total =ls_pair_total + ls_self_total + ls_surface_total;
  
  nonbonded_total = lj_total + crf_total + self_total + 
          ls_realspace_total + ls_kspace_total + ls_self_total + ls_surface_total +
          ls_a_term_total;
  bonded_total = bond_total + angle_total + dihedral_total + improper_total
                 + crossdihedral_total;
  potential_total = nonbonded_total + bonded_total;
  
  special_total = posrest_total + distanceres_total + dihrest_total
    + constraints_total + jvalue_total + xray_total + external_total 
    + eds_vr + leus_total + sasa_total + sasa_volume_total;

  total = potential_total + kinetic_total + special_total;

/*
#ifdef HAVE_ISNAN
  if (std::isnan(total)){
    io::messages.add("total energy is NaN", "energy", io::message::error);
    return E_NAN;
   }
#endif
*/
  if (math::gisnan(total)) {
    io::messages.add("total energy is NaN", "energy", io::message::error);
    return E_NAN;    
  }

  return 0;
}

double configuration::Energy::get_energy_by_index(const unsigned int & index) {
  switch(index) {
    case 1 : return total;
    case 2 : return potential_total;
    case 3 : return bond_total;
    case 4 : return angle_total;
    case 5 : return improper_total;
    case 6 : return dihedral_total;
    case 7 : return crossdihedral_total;
    case 8 : return bonded_total;
    case 9 : return nonbonded_total;
    case 10 : return lj_total;
    case 11 : return crf_total;
    case 12 : return ls_total;
    case 13 : return ls_pair_total;
    case 14 : return ls_realspace_total;
    case 15 : return ls_kspace_total;
    case 16 : return ls_self_total;
    case 17 : return ls_surface_total;
    case 18 : return ls_a_term_total;
    case 19 : return special_total;
    case 20 : return posrest_total;
    case 21 : return distanceres_total;
    case 22 : return dihrest_total;
    case 23 : return jvalue_total;
    case 24 : return xray_total;
    case 25 : return leus_total;
    case 26 : return constraints_total;
    case 27 : return entropy_term;
    case 28 : return external_total;
    case 29 : return self_total;
    case 30 : return eds_vr;
    case 31 : return sasa_total;
    case 32 : return sasa_volume_total;
  }
  return 0.0;
}
