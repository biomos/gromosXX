/**
 * @file solute.tcc
 * inline methods for solute-topology.
 */


/**
 * atom information accessor.
 */
inline 
simulation::Solute::atom_struct & 
simulation::Solute::atom(size_t i)
{
  assert(i < m_atom.size());
  return m_atom[i];
}

/**
 * all atom information accessor.
 */
inline 
std::vector<simulation::Solute::atom_struct> & 
simulation::Solute::atoms()
{
  return m_atom;
}

/**
 * add a solute atom.
 */
inline void 
simulation::Solute::add_atom(std::string name, int residue_nr)
{
  atom_struct s;
  s.name = name;
  s.residue_nr = residue_nr;
  
  m_atom.push_back(s);
  ++m_num_atoms;
}

/**
 * bonds accessor.
 */
inline simulation::Bond &
simulation::Solute::bonds()
{
  return m_bond;
}

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
inline void 
simulation::Solute::add_bond_length_constraints(std::vector<interaction::bond_type_struct> const &param)
{
  Bond bonds;
  Bond::iterator it = m_bond.begin();
  for( ; it.neol(); ++it){
    add_distance_constraint(it.i(), it.j(), param[it.type()].r0);
  }
  m_bond = bonds;
}
    
/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
inline void
simulation::Solute
::add_bond_length_constraints(int iac,
			     std::vector<int> const &atom_iac,
			     std::vector<interaction::bond_type_struct> const &param)
{
  Bond bonds;
  Bond::iterator it = m_bond.begin();
  for( ; it.neol(); ++it){
    if(atom_iac[it.i()] == iac || atom_iac[it.j()] == iac)
      add_distance_constraint(it.i(), it.j(), param[it.type()].r0);
    else
      bonds.add(it.i(), it.j(), it.type());
  }
  m_bond = bonds;
}
    
/**
 * add bonds connecting an atom of mass mass to the
 * constraint vector and remove from the bond vector.
 */
inline void
simulation::Solute::add_bond_length_constraints(double mass,
					       math::SArray const &atom_mass,
					       std::vector<interaction::bond_type_struct> const &param)
{
  Bond bonds;
  Bond::iterator it = m_bond.begin();
  for( ; it.neol(); ++it){
    if(atom_mass(it.i()) == mass || atom_mass(it.j()) == mass)
      add_distance_constraint(it.i(), it.j(), param[it.type()].r0);
    else
      bonds.add(it.i(), it.j(), it.type());
  }
  m_bond = bonds;
}


/**
 * angles accessor.
 */
inline simulation::Angle &
simulation::Solute::angles()
{
  return m_angle;
}

/**
 * improper_dihedrals accessor.
 */
inline simulation::Improper_dihedral &
simulation::Solute::improper_dihedrals()
{
  return m_improper_dihedral;
}

/**
 * dihedrals accessor.
 */
inline simulation::Dihedral &
simulation::Solute::dihedrals()
{
  return m_dihedral;
}
