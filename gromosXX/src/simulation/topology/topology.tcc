/**
 * @file topology.tcc
 * inline methods definition
 */

/**
 * Constructor
 */
inline simulation::topology::topology()
  : m_num_solute_atoms(0), 
    m_mass(0),
    m_charge(0)
{
}

/**
 * mass accessor
 */
math::SArray & simulation::topology::mass()
{
  return m_mass;
}

/**
 * const mass accessor
 */
math::SArray const & simulation::topology::mass()const
{
  return m_mass;
}

/**
 * charge accessor
 */
math::SArray & simulation::topology::charge()
{
  return m_charge;
}

/**
 * const charge accessor
 */
math::SArray const & simulation::topology::charge()const
{
  return m_charge;
}

/**
 * the number of solute atoms
 */
inline size_t simulation::topology::num_solute_atoms()const
{
  return m_num_solute_atoms;
}

/**
 * set the capacity of solute atoms byresizeing
 * the apropriate arrays.
 */
inline void simulation::topology::solute_atoms_capacity(size_t atoms)
{
  m_mass.resizeAndPreserve(atoms);
  m_charge.resizeAndPreserve(atoms);
  m_exclusion.resize(atoms);
  m_one_four_pair.resize(atoms);
  m_all_exclusion.resize(atoms);
}

/**
 * set the number of solute atoms
 */
inline void simulation::topology::num_solute_atoms(size_t atoms)
{
  assert(unsigned(m_mass.size()) <= atoms);
  m_num_solute_atoms = atoms;
}

/**
 * bond accessor
 */
inline simulation::bond & simulation::topology::bonds()
{
  return m_bonds;
}

namespace simulation
{
  /**
   * output information about the topology.
   */
  inline std::ostream & operator<<(std::ostream &os, topology &topo)
  {
    os << "a topology";
    return os;
  }
}

/**
 * add a solute atom to the topology.
 * if the arrays are too small they will be increased.
 * if adding multiple solute atoms, first call solute_atoms_capacity...
 */
inline void simulation::topology::add_solute_atom(std::string name, int residue_nr,
						  int iac, double mass,
						  double charge, bool chargegroup,
						  std::set<int> exclusions,
						  std::set<int> one_four_pairs)
{

  if (unsigned(m_mass.size()) < m_num_solute_atoms + 1){
    solute_atoms_capacity(m_num_solute_atoms+1);
  }
  
  soluteatom().add(name, residue_nr, iac);

  topology::mass()(m_num_solute_atoms) = mass;
  topology::charge()(m_num_solute_atoms) = charge;

  if (chargegroup) m_chargegroup.push_back(m_num_solute_atoms+1);
  
  m_exclusion[m_num_solute_atoms] = exclusions;
  m_one_four_pair[m_num_solute_atoms] = one_four_pairs;
  
  set_union(exclusions.begin(), exclusions.end(),
	    one_four_pairs.begin(), one_four_pairs.end(),
	    inserter(m_all_exclusion[m_num_solute_atoms], 
		     m_all_exclusion[m_num_solute_atoms].end())
	    );

  ++m_num_solute_atoms;  

}
