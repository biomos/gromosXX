/**
 * @file solvent.tcc
 * inline methods for solvent topology.
 */

/**
 * add a solvent atom.
 */
inline void simulation::Solvent
::add_atom(std::string name, int res_nr, int iac,
	   double mass, double charge)
{
  atom_struct s;
  s.name = name;
  s.residue_nr = res_nr;
  s.iac = iac;
  s.mass = mass;
  s.charge = charge;

  m_atom.push_back(s);
  ++m_num_atoms;
  
}

/**
 * accessor - atom i.
 */
inline simulation::Solvent::atom_struct &
simulation::Solvent::atom(size_t i)
{
  assert(i < m_atom.size());
  return m_atom[i];
}

/**
 * accessor - vector.
 */
inline
std::vector<simulation::Solvent::atom_struct> &
simulation::Solvent::atoms()
{
  return m_atom;
}
