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
inline simulation::bond &
simulation::Solute::bonds()
{
  return m_bond;
}

/**
 * add all bonds to the solute constraint vector and
 * remove them from the bond vector.
 */
inline void 
simulation::Solute::add_bond_length_constraints()
{
}
    
/**
 * add bonds connecting an atom of type iac to the
 * constraint vector and remove from the bond vector.
 */
inline void
simulation::Solute::add_bond_length_constraint(int iac)
{
}
    
/**
 * add bonds connecting an atom of mass mass to the
 * constraint vector and remove from the bond vector.
 */
inline void
simulation::Solute::add_bond_length_constraint(double mass)
{
}


