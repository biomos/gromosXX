/**
 * @file topology.tcc
 * inline methods definition
 */

/**
 * Constructor
 */
inline simulation::topology::topology()
  : m_num_solute_atoms(0),
    m_num_solvent_atoms(0),
    m_mass(0),
    m_charge(0)
{
}

/**
 * integer atom code accessor.
 */
inline int simulation::topology::iac(int i)
{
  return m_iac[i];
}

/**
 * mass accessor
 */
inline math::SArray & simulation::topology::mass()
{
  return m_mass;
}

/**
 * const mass accessor
 */
inline math::SArray const & simulation::topology::mass()const
{
  return m_mass;
}

/**
 * charge accessor
 */
inline math::SArray & simulation::topology::charge()
{
  return m_charge;
}

/**
 * const charge accessor
 */
inline math::SArray const & simulation::topology::charge()const
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
inline void simulation::topology::resize(size_t atoms)
{
  // if you want to shrink, first change num_atoms...
  assert(atoms >= num_solute_atoms() + num_solvent_atoms());

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
  assert(num_solvent_atoms() == 0);
  assert(unsigned(m_mass.size()) == atoms);

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
    resize(m_num_solute_atoms+1);
  }
  
  soluteatoms().add(name, residue_nr);

  topology::mass()(m_num_solute_atoms) = mass;
  topology::charge()(m_num_solute_atoms) = charge;

  if (chargegroup) m_chargegroup.push_back(m_num_solute_atoms+1);
  
  m_iac.push_back(iac);

  m_exclusion[m_num_solute_atoms] = exclusions;
  m_one_four_pair[m_num_solute_atoms] = one_four_pairs;
  
  set_union(exclusions.begin(), exclusions.end(),
	    one_four_pairs.begin(), one_four_pairs.end(),
	    inserter(m_all_exclusion[m_num_solute_atoms], 
		     m_all_exclusion[m_num_solute_atoms].end())
	    );

  ++m_num_solute_atoms;  

}

/**
 * soluteatom accessor.
 */
inline simulation::soluteatom & simulation::topology::soluteatoms()
{
  return m_soluteatoms;
}

/**
 * solvent accessor.
 */
inline simulation::solvent & simulation::topology::solvents(size_t i)
{
  assert(i < m_solvents.size());
  return m_solvents[i];
}

/**
 * number of solvents.
 */
inline size_t simulation::topology::num_solvents()const
{
  return m_num_solvent_molecules.size();
}

/**
 * add a solvent.
 */
inline void simulation::topology::add_solvent(solvent solv)
{
  m_solvents.push_back(solv);
}

/**
 * add solvent to the simulation.
 */
inline void simulation::topology::solvate(size_t solv, size_t num_molecules)
{
  // only add in the correct order!
  assert(solv == m_num_solvent_atoms.size());
  assert(solv < m_solvents.size());

  int n = num_solute_atoms() + num_solvent_atoms();

  m_num_solvent_molecules.push_back(num_molecules);
  m_num_solvent_atoms.push_back(num_molecules * m_solvents[solv].num_atoms());
  
  resize(num_solute_atoms() + num_solvent_atoms());

  // add to iac, mass, charge
  for(size_t i=0; i<num_molecules; ++i){
    for(size_t j=0; j<m_solvents[solv].num_atoms(); ++j, ++n){
      m_iac.push_back(m_solvents[solv].atom(j).iac);
      m_mass(n) = m_solvents[solv].atom(j).mass;
      m_charge(n) = m_solvents[solv].atom(j).charge;
      // no exclusions or 1-4 interactions for solvent ?!
    }
  }
  
}

/**
 * number of solvent molecules.
 */
inline size_t simulation::topology::num_solvent_molecules(size_t i)const
{
  assert(i < m_num_solvent_molecules.size());
  return m_num_solvent_molecules[i];
}

/**
 * total solvent atoms accessor
 */
inline size_t simulation::topology::num_solvent_atoms()const
{
  size_t n = 0;
  for(std::vector<size_t>::const_iterator it = m_num_solvent_atoms.begin(),
	to = m_num_solvent_atoms.end();
      it != to; ++it)
    n += *it;
  return n;
}

/**
 * solvent atoms of solvent i
 */
inline size_t simulation::topology::num_solvent_atoms(size_t i)const
{
  assert(i<m_num_solvent_atoms.size());
  return m_num_solvent_atoms[i];
}

/**
 * residue name accessor.
 */
inline std::vector<std::string> & simulation::topology::residue_name()
{
  return m_residue_name;
}
