/**
 * @file soluteatom.tcc
 * inline methods for soluteatom-topology.
 */

/**
 * add a solute atom.
 */
inline void simulation::soluteatom::add(std::string name, int residue_nr, int iac)
{
  soluteatom_struct s;
  s.name = name;
  s.residue_nr = residue_nr;
  s.iac = iac;
  
  m_information.push_back(s);
}

/**
 * accessor.
 */
inline 
simulation::soluteatom::soluteatom_struct & 
simulation::soluteatom::operator()(size_t i)
{
  return m_information[i];
}
