/**
 * @file solute.tcc
 * inline methods for solute-topology.
 */

/**
 * Constructor
 */
inline
simulation::Solute::Solute()
  : compound()
{
}

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
inline std::vector<simulation::Bond> &
simulation::Solute::bonds()
{
  return m_bond;
}


/**
 * angles accessor.
 */
inline std::vector<simulation::Angle> &
simulation::Solute::angles()
{
  return m_angle;
}

/**
 * improper_dihedrals accessor.
 */
inline std::vector<simulation::Improper_Dihedral> &
simulation::Solute::improper_dihedrals()
{
  return m_improper_dihedral;
}

/**
 * dihedrals accessor.
 */
inline std::vector<simulation::Dihedral> &
simulation::Solute::dihedrals()
{
  return m_dihedral;
}
