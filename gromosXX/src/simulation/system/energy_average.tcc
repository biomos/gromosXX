/**
 * @file energy_average.tcc
 * implements member methods.
 */

inline
simulation::Energy_Average::Energy_Average()
  : m_average(),
    m_square_average(),
    m_time(0)
{
}

inline void
simulation::Energy_Average::zero()
{
  m_average.zero();
  m_square_average.zero();
  m_time = 0.0;
}

inline void
simulation::Energy_Average::resize(size_t const energy_groups, size_t const multi_baths)
{
  m_average.resize(energy_groups, multi_baths);
  m_square_average.resize(energy_groups, multi_baths);
  m_time = 0.0;
}

inline void
simulation::Energy_Average::update(simulation::Energy const &e, double const dt)
{
  // totals
  m_average.total += dt * e.total;
  
  m_square_average.total += dt * e.total * e.total;
  
  m_average.kinetic_total += dt * e.kinetic_total;
  m_square_average.kinetic_total += dt * e.kinetic_total * e.kinetic_total;
  
  m_average.potential_total += dt * e.potential_total;
  m_square_average.potential_total += dt * e.potential_total * e.potential_total;
  
  // physical interaction totals
  m_average.bond_total += dt * e.bond_total;
  m_square_average.bond_total = dt * e.bond_total * e.bond_total;
  
  m_average.angle_total += dt * e.angle_total;
  m_square_average.angle_total += dt * e.angle_total * e.angle_total;
  
  m_average.improper_total += dt * e.improper_total;
  m_square_average.improper_total += dt * e.improper_total * e.improper_total;
  
  m_average.dihedral_total += dt * e.dihedral_total;
  m_square_average.dihedral_total += dt * e.dihedral_total * e.dihedral_total;
  
  m_average.bonded_total += dt * e.bonded_total;
  m_square_average.bonded_total += dt * e.bonded_total * e.bonded_total;
  
  m_average.nonbonded_total += dt * e.nonbonded_total;
  m_square_average.nonbonded_total += dt * e.nonbonded_total * e.nonbonded_total;
  
  m_average.lj_total += dt * e.lj_total;
  m_square_average.lj_total += dt * e.lj_total * e.lj_total;
  
  m_average.crf_total += dt * e.crf_total;
  m_square_average.crf_total += dt * e.crf_total * e.crf_total;
  
  m_average.special_total += dt * e.special_total;
  m_square_average.special_total += dt * e.special_total * e.special_total;
  
  // kinetic energies of the baths
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){
    m_average.kinetic_energy[i] += dt * e.kinetic_energy[i];
    m_square_average.kinetic_energy[i] += dt * e.kinetic_energy[i] * e.kinetic_energy[i];
  }
  
  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){
    m_average.bond_energy[i] += dt * e.bond_energy[i];
    m_square_average.bond_energy[i] = dt * e.bond_energy[i] * e.bond_energy[i];
  
    m_average.angle_energy[i] += dt * e.angle_energy[i];
    m_square_average.angle_energy[i] += dt * e.angle_energy[i] * e.angle_energy[i];
  
    m_average.improper_energy[i] += dt * e.improper_energy[i];
    m_square_average.improper_energy[i] += dt * e.improper_energy[i] * e.improper_energy[i];
  
    m_average.dihedral_energy[i] += dt * e.dihedral_energy[i];
    m_square_average.dihedral_energy[i] += dt * e.dihedral_energy[i] * e.dihedral_energy[i];

    // and the nonbonded groups
    for(size_t j=0; j < e.lj_energy.size(); ++j){
      m_average.lj_energy[i][j] += dt * e.lj_energy[i][j];
      m_square_average.lj_energy[i][j] += dt * e.lj_energy[i][j] * e.lj_energy[i][j];
  
      m_average.crf_energy[i][j] += dt * e.crf_energy[i][j];
      m_square_average.crf_energy[i][j] += dt * e.crf_energy[i][j] * e.crf_energy[i][j];
      
    }
  
  }

  m_time += dt;
  
}

inline void
simulation::Energy_Average::average(simulation::Energy &energy, simulation::Energy &fluctuation)
{
  // only at the end of a run, so create a temporary
  Energy &e = energy;
  Energy &f = fluctuation;

  e.resize(m_average.bond_energy.size(), m_average.kinetic_energy.size());
  f.resize(m_average.bond_energy.size(), m_average.kinetic_energy.size());

  // totals
  e.total = m_average.total / m_time;
  f.total = sqrt((m_square_average.total - m_average.total * m_average.total / m_time) / m_time);

  e.kinetic_total = m_average.kinetic_total / m_time;
  f.kinetic_total = sqrt((m_square_average.kinetic_total -
			  m_average.kinetic_total * m_average.kinetic_total / m_time) / m_time);

  e.potential_total = m_average.potential_total / m_time;
  f.potential_total = sqrt((m_square_average.potential_total - 
			    m_average.potential_total * m_average.potential_total / m_time) / m_time);
  
  // physical interaction totals
  e.bond_total = m_average.bond_total / m_time;
  f.bond_total = sqrt((m_square_average.bond_total - 
		       m_average.bond_total * m_average.bond_total / m_time) / m_time);

  e.angle_total = m_average.angle_total / m_time;
  f.angle_total = sqrt((m_square_average.angle_total -
			m_average.angle_total * m_average.angle_total / m_time) / m_time);

  e.improper_total = m_average.improper_total / m_time;
  f.improper_total = sqrt((m_square_average.improper_total - 
			   m_average.improper_total * m_average.improper_total / m_time) / m_time);

  e.dihedral_total = m_average.dihedral_total / m_time;
  f.dihedral_total = sqrt((m_square_average.dihedral_total -
			   m_average.dihedral_total * m_average.dihedral_total / m_time) / m_time);

  e.bonded_total = m_average.bonded_total / m_time;
  f.bonded_total = sqrt((m_square_average.bonded_total -
			 m_average.bonded_total * m_average.bonded_total / m_time) / m_time);

  e.nonbonded_total = m_average.nonbonded_total / m_time;
  f.nonbonded_total = sqrt((m_square_average.nonbonded_total -
			    m_average.nonbonded_total * m_average.nonbonded_total / m_time) / m_time);
  
  e.lj_total = m_average.lj_total / m_time;
  f.lj_total = sqrt((m_square_average.lj_total -
		     m_average.lj_total * m_average.lj_total / m_time) / m_time);

  e.crf_total = m_average.crf_total / m_time;
  f.crf_total = sqrt((m_square_average.crf_total -
		      m_average.crf_total * m_average.crf_total / m_time) / m_time);

  e.special_total = m_average.special_total / m_time;
  f.special_total = sqrt((m_square_average.special_total -
			  m_average.special_total * m_average.special_total / m_time) / m_time);
  
  // kinetic energies of the baths
  for(size_t i=0; i < e.kinetic_energy.size(); ++i){
    e.kinetic_energy[i] = m_average.kinetic_energy[i] / m_time;
    f.kinetic_energy[i] = sqrt((m_square_average.kinetic_energy[i] -
				m_average.kinetic_energy[i] * m_average.kinetic_energy[i] / m_time)
			       / m_time);
  }
  
  // the energy groups
  for(size_t i=0; i < e.bond_energy.size(); ++i){
    e.bond_energy[i] = m_average.bond_energy[i] / m_time;
    f.bond_energy[i] = sqrt((m_square_average.bond_energy[i] -
			     m_average.bond_energy[i] * m_average.bond_energy[i] / m_time) / m_time);

    e.angle_energy[i] = m_average.angle_energy[i] / m_time;
    f.angle_energy[i] = sqrt((m_square_average.angle_energy[i] -
			 m_average.angle_energy[i] * m_average.angle_energy[i] / m_time) / m_time);

    e.improper_energy[i] = m_average.improper_energy[i] / m_time;
    f.improper_energy[i] = sqrt((m_square_average.improper_energy[i]
				 - m_average.improper_energy[i] * m_average.improper_energy[i] / m_time)
				/ m_time);

    e.dihedral_energy[i] = m_average.dihedral_energy[i] / m_time;
    f.dihedral_energy[i] = sqrt((m_square_average.dihedral_energy[i] -
				 m_average.dihedral_energy[i] * m_average.dihedral_energy[i] / m_time)
				/ m_time);

    // and the nonbonded groups
    for(size_t j=0; j < e.lj_energy.size(); ++j){
      e.lj_energy[i][j] = m_average.lj_energy[i][j] / m_time;
      f.lj_energy[i][j] = sqrt((m_square_average.lj_energy[i][j] -
				m_average.lj_energy[i][j] * m_average.lj_energy[i][j] / m_time) / m_time);

      e.crf_energy[i][j] = m_average.crf_energy[i][j] / m_time;
      f.crf_energy[i][j] = sqrt((m_square_average.crf_energy[i][j] -
				 m_average.crf_energy[i][j] * m_average.crf_energy[i][j] / m_time)
				/ m_time);
    }
  
  }

}
