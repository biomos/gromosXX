/**
 * @file energy.tcc
 * implements the energy methods.
 */

inline void simulation::Energy::zero()
{
  bond_energy.assign(bond_energy.size(), 0.0);
  angle_energy.assign(angle_energy.size(), 0.0);
  improper_energy.assign(improper_energy.size(), 0.0);
  dihedral_energy.assign(dihedral_energy.size(), 0.0);
  
  lj_energy.assign(lj_energy.size(), 
		   std::vector<double>(lj_energy.size(), 0.0));
  crf_energy.assign(crf_energy.size(), 
		    std::vector<double>(crf_energy.size(), 0.0));
}

inline void simulation::Energy::resize(size_t s)
{
  bond_energy.resize(s);
  angle_energy.resize(s);
  improper_energy.resize(s);
  dihedral_energy.resize(s);
  
  lj_energy.resize(s);
  crf_energy.resize(s);

  for(size_t i=0; i<s; ++i){
    lj_energy[i].resize(s);
    crf_energy[i].resize(s);
  }

  zero();
  
}

