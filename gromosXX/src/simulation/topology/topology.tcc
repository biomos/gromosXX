/**
 * @file topology.tcc
 * inline methods definition
 */

/**
 * Constructor
 */
inline simulation::topology::topology()
  : m_mass(0),
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
 * solute accessor.
 */
inline simulation::Solute &
simulation::topology::solute()
{
  return m_solute;
}

/**
 * const solute accessor.
 */
inline simulation::Solute const &
simulation::topology::solute()const
{
  return m_solute;
}

/**
 * the number of solute atoms
 */
inline size_t simulation::topology::num_solute_atoms()const
{
  return solute().num_atoms();
}

/**
 * set the capacity of solute atoms by resizeing
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
 * add a solute atom to the topology.
 * if the arrays are too small they will be increased.
 * if adding multiple solute atoms, first call solute_atoms_capacity...
 */
inline void simulation::topology
::add_solute_atom(std::string name, int residue_nr,
		  int iac, double mass,
		  double charge, bool chargegroup,
		  std::set<int> exclusions,
		  std::set<int> one_four_pairs)
{

  if (unsigned(m_mass.size()) < num_solute_atoms() + 1){
    resize(num_solute_atoms()+1);
  }
  
  topology::mass()(num_solute_atoms()) = mass;
  topology::charge()(num_solute_atoms()) = charge;

  if (chargegroup) m_chargegroup.push_back(num_solute_atoms()+1);
  
  m_iac.push_back(iac);

  m_exclusion[num_solute_atoms()] = exclusions;
  m_one_four_pair[num_solute_atoms()] = one_four_pairs;
  
  std::set_union(exclusions.begin(), exclusions.end(),
		 one_four_pairs.begin(), one_four_pairs.end(),
		 std::inserter(m_all_exclusion[num_solute_atoms()], 
			       m_all_exclusion[num_solute_atoms()].end())
		 );

  // this increases num_solute_atoms()
  solute().add_atom(name, residue_nr);

}

/**
 * solvent accessor.
 */
inline simulation::Solvent & simulation::topology::solvent(size_t i)
{
  assert(i < m_solvent.size());
  return m_solvent[i];
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
inline void simulation::topology::add_solvent(Solvent solv)
{
  m_solvent.push_back(solv);
}

/**
 * add solvent molecules to the simulation (system).
 */
inline void simulation::topology::solvate(size_t solv, size_t num_molecules)
{
  // only add in the correct order!
  assert(solv == m_num_solvent_atoms.size());
  assert(solv < m_solvent.size());

  int n = num_solute_atoms() + num_solvent_atoms();

  m_num_solvent_molecules.push_back(num_molecules);
  m_num_solvent_atoms.push_back(num_molecules * m_solvent[solv].num_atoms());
  
  resize(num_solute_atoms() + num_solvent_atoms());

  // add to iac, mass, charge
  for(size_t i=0; i<num_molecules; ++i){
    for(size_t j=0; j<m_solvent[solv].num_atoms(); ++j, ++n){

      m_iac.push_back(m_solvent[solv].atom(j).iac);
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;
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
 * total number of solvent atoms.
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
 * solvent atoms of solvent i (*molecules).
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

/**
 * all exclusions for atom i.
 */
inline std::set<int> & simulation::topology::all_exclusion(size_t i)
{
  return m_all_exclusion[i];
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

