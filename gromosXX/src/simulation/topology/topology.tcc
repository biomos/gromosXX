/**
 * @file topology.tcc
 * inline methods definition
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE topology

#include "../../debug.h"

/**
 * Constructor
 */
inline simulation::Topology::Topology()
  : m_mass(0),
    m_charge(0),
    m_num_solute_chargegroups(0)
{
  m_chargegroup.push_back(0);
  m_molecule.push_back(0);
}

/**
 * integer atom code accessor.
 */
inline int simulation::Topology::iac(size_t const i)const
{
  assert(i < m_iac.size());
  return m_iac[i];
}

/**
 * mass accessor
 */
inline math::SArray & simulation::Topology::mass()
{
  return m_mass;
}

/**
 * const mass accessor
 */
inline math::SArray const & simulation::Topology::mass()const
{
  return m_mass;
}

/**
 * charge accessor
 */
inline math::SArray & simulation::Topology::charge()
{
  return m_charge;
}

/**
 * const charge accessor
 */
inline math::SArray const & simulation::Topology::charge()const
{
  return m_charge;
}

/**
 * solute accessor.
 */
inline simulation::Solute &
simulation::Topology::solute()
{
  return m_solute;
}

/**
 * const solute accessor.
 */
inline simulation::Solute const &
simulation::Topology::solute()const
{
  return m_solute;
}

/**
 * the number of atoms
 */
inline size_t simulation::Topology::num_atoms()const
{
  return num_solute_atoms() + num_solvent_atoms();
}

/**
 * the number of solute atoms
 */
inline size_t simulation::Topology::num_solute_atoms()const
{
  return solute().num_atoms();
}

/**
 * set the capacity of solute atoms by resizeing
 * the apropriate arrays.
 */
inline void simulation::Topology::resize(size_t const atoms)
{
  // if you want to shrink, first change num_atoms...
  assert(atoms >= num_solute_atoms() + num_solvent_atoms());

  m_mass.resizeAndPreserve(atoms);
  m_charge.resizeAndPreserve(atoms);
  m_exclusion.resize(atoms);
  m_one_four_pair.resize(atoms);
  m_all_exclusion.resize(atoms);
  
  m_iac.resize(atoms);
  // chargegroups???
}

/**
 * add a solute atom to the topology.
 * if the arrays are too small they will be increased.
 * if adding multiple solute atoms, first call solute_atoms_capacity...
 */
inline void simulation::Topology
::add_solute_atom(std::string name, int residue_nr,
		  int iac, double mass,
		  double charge, bool chargegroup,
		  std::set<int> exclusions,
		  std::set<int> one_four_pairs)
{

  if (unsigned(m_mass.size()) < num_solute_atoms() + 1){
    resize(num_solute_atoms()+1);
  }
  
  Topology::mass()(num_solute_atoms()) = mass;
  Topology::charge()(num_solute_atoms()) = charge;

  if (chargegroup){
    m_chargegroup.push_back(num_solute_atoms()+1);
    ++m_num_solute_chargegroups;
  }
  
  DEBUG(15, "iac[" << num_solute_atoms() << "] = " << iac);

  m_iac[num_solute_atoms()] = iac;

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
inline simulation::Solvent & simulation::Topology::solvent(size_t i)
{
  assert(i < m_solvent.size());
  return m_solvent[i];
}
/**
 * const solvent accessor.
 */
inline simulation::Solvent const & simulation::Topology::solvent(size_t i)const
{
  assert(i < m_solvent.size());
  return m_solvent[i];
}

/**
 * number of solvents.
 */
inline size_t simulation::Topology::num_solvents()const
{
  return m_num_solvent_molecules.size();
}

/**
 * add a solvent.
 */
inline void simulation::Topology::add_solvent(Solvent solv)
{
  m_solvent.push_back(solv);
}

/**
 * add solvent molecules to the simulation (system).
 */
inline void simulation::Topology::solvate(size_t solv, size_t num_molecules)
{
  // only add in the correct order!
  assert(solv == m_num_solvent_atoms.size());
  assert(solv < m_solvent.size());

  int n = num_solute_atoms() + num_solvent_atoms();

  m_num_solvent_molecules.push_back(num_molecules);
  m_num_solvent_atoms.push_back(num_molecules * m_solvent[solv].num_atoms());
  
  DEBUG(5, "solvate: solvent atoms: " << num_solvent_atoms());
  DEBUG(10, "solvate: total atoms: " << num_solute_atoms() + num_solvent_atoms());
  
  resize(num_solute_atoms() + num_solvent_atoms());

  // add to iac, mass, charge
  for(size_t i=0; i<num_molecules; ++i){
    for(size_t j=0; j<m_solvent[solv].num_atoms(); ++j, ++n){

      DEBUG(15, "iac[" << n << "]=" << m_solvent[solv].atom(j).iac);
      DEBUG(15, "charge[" << n << "]=" << m_solvent[solv].atom(j).charge);
      
      m_iac[n] = m_solvent[solv].atom(j).iac;
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;
      // no exclusions or 1-4 interactions for solvent ?!

    }

    // add to the chargegroups
    DEBUG(11, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    m_molecule.push_back(n);

  }
    
}

/**
 * number of solvent molecules.
 */
inline size_t simulation::Topology::num_solvent_molecules(size_t i)const
{
  assert(i < m_num_solvent_molecules.size());
  return m_num_solvent_molecules[i];
}

/**
 * total number of solvent atoms.
 */
inline size_t simulation::Topology::num_solvent_atoms()const
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
inline size_t simulation::Topology::num_solvent_atoms(size_t i)const
{
  assert(i<m_num_solvent_atoms.size());
  return m_num_solvent_atoms[i];
}

/**
 * residue name accessor.
 */
inline std::vector<std::string> & simulation::Topology::residue_name()
{
  return m_residue_name;
}

/**
 * all exclusions for atom i. Exclusions and 1,4 interactions.
 */
inline std::set<int> & simulation::Topology::all_exclusion(size_t const i)
{
  assert(i < m_all_exclusion.size());
  return m_all_exclusion[i];
}

inline std::set<int> const & simulation::Topology::all_exclusion(size_t const i)const
{
  assert(i < m_all_exclusion.size());
  return m_all_exclusion[i];
}

/**
 * exclusions for atom i.
 */
inline std::set<int> & simulation::Topology::exclusion(size_t const i)
{
  assert(i < m_exclusion.size());
  return m_exclusion[i];
}
/**
 * all exclusions
 */
inline std::vector<std::set<int> > & simulation::Topology::exclusion()
{
  return m_exclusion;
}


/**
 * all one four pairs for atom i.
 */
inline std::set<int> & simulation::Topology::one_four_pair(size_t const i)
{
  assert(i < m_one_four_pair.size());
  return m_one_four_pair[i];
}
/**
 * all one four pairs.
 */
inline std::vector<std::set<int> > & simulation::Topology::one_four_pair()
{
  return m_one_four_pair;
}

/**
 * iterator over the chargegrops
 */
inline simulation::chargegroup_iterator 
simulation::Topology::chargegroup_begin()const
{
  return chargegroup_iterator(m_chargegroup.begin());
}

/**
 * end of the chargegroup iterator.
 */
inline simulation::chargegroup_iterator
simulation::Topology::chargegroup_end()const
{
  return chargegroup_iterator(m_chargegroup.end()-1);
}

/**
 * iterator over the molecules
 */
inline simulation::Molecule_Iterator 
simulation::Topology::molecule_begin()
{
  return Molecule_Iterator(m_molecule.begin());
}

/**
 * end of the molecule iterator.
 */
inline simulation::Molecule_Iterator
simulation::Topology::molecule_end()
{
  return Molecule_Iterator(m_molecule.end()-1);
}

/**
 * the number of chargegroups present.
 */
inline size_t simulation::Topology::num_chargegroups()const
{
  return m_chargegroup.size()-1;
}

/**
 * the number of solute chargegroups.
 */
inline size_t simulation::Topology::num_solute_chargegroups()const
{
  return m_num_solute_chargegroups;
}

/**
 * the molecule indices.
 */
inline std::vector<size_t> & simulation::Topology::molecules()
{
  return m_molecule;
}

/**
 * const energy group accessor.
 */
inline std::vector<size_t> const & simulation::Topology::energy_groups()const
{
  return m_energy_group;
}

/**
 * energy group accessor.
 */
inline std::vector<size_t> & simulation::Topology::energy_groups()
{
  return m_energy_group;
}

/**
 * const atom energy group accessor.
 */
inline std::vector<size_t> const & 
simulation::Topology::atom_energy_group()const
{
  return m_atom_energy_group;
}

/**
 * atom energy group accessor.
 */
inline std::vector<size_t> & 
simulation::Topology::atom_energy_group()
{
  return m_atom_energy_group;
}

/**
 * atom energy group accessor.
 */
inline const size_t 
simulation::Topology::atom_energy_group(size_t i)const
{
  assert(i < m_atom_energy_group.size());
  return m_atom_energy_group[i];
}

/**
 * calculate constraint degrees of freedom.
 */
inline void
simulation::Topology::
calculate_constraint_dof(simulation::Multibath &multibath)const
{
  // substract constraints
  std::vector<compound::distance_constraint_struct>::const_iterator 
    c_it = solute().distance_constraints().begin(),
    c_to = solute().distance_constraints().end();
  
  size_t com_bath_i, ir_bath_i, com_bath_j, ir_bath_j;

  for( ; c_it != c_to; ++c_it){
    
    DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
    multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
    multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);

    multibath[ir_bath_i].dof -= 0.5;
    multibath[ir_bath_j].dof -= 0.5;

    multibath[ir_bath_i].ir_dof -= 0.5;
    multibath[ir_bath_j].ir_dof -= 0.5;

    multibath[ir_bath_i].solute_constr_dof += 0.5;
    multibath[ir_bath_j].solute_constr_dof += 0.5;

  }
  
  for(size_t i=0; i<multibath.size(); ++i){
    DEBUG(7, "dof           " << multibath[i].dof);
    DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
  }

  // solvent constraints
  int index = num_solute_atoms();
  for(size_t s=0; s < num_solvents(); ++s){
    
    for(size_t m=0; m < num_solvent_molecules(s); ++m){

      c_it = solvent(s).distance_constraints().begin();
      c_to = solvent(s).distance_constraints().end();
      
      for( ; c_it != c_to; ++c_it){

	multibath.in_bath(c_it->i + index, com_bath_i, ir_bath_i);
	multibath.in_bath(c_it->j + index, com_bath_j, ir_bath_j);
	
	multibath[ir_bath_i].dof -= 0.5;
	multibath[ir_bath_j].dof -= 0.5;

	multibath[ir_bath_i].ir_dof -= 0.5;
	multibath[ir_bath_j].ir_dof -= 0.5;
	
	multibath[ir_bath_i].solvent_constr_dof += 0.5;
	multibath[ir_bath_j].solvent_constr_dof += 0.5;
	
      }
      
      index += solvent(s).num_atoms();
      
    }
    
  }

}

/**
 * check state
 */
inline int 
simulation::Topology::check_state()const
{
  int result = 0;

  // check that we have an energy group for every atom
  if (m_atom_energy_group.size() != num_atoms()){
    io::messages.add("not every atom has an energy group index",
		     "Topology::check_state", io::message::error);
    ++result;
  }
  for(std::vector<size_t>::const_iterator it = m_atom_energy_group.begin(),
	to = m_atom_energy_group.end(); it != to; ++it){
    if (*it >= m_energy_group.size()){
      io::messages.add("energy group index of atom too large",
		       "Topology::check_state", io::message::error);
      ++result;
    }
  }

  return result;
}

/**
 * return a lambda for a nonperturbed system
 */
inline const double simulation::Topology::lambda()const
{
  return 0.0;
}

inline const int simulation::Topology::nlam()const
{
  return 1;
}

namespace simulation
{
  /**
   * output information about the topology.
   */
  inline std::ostream & operator<<(std::ostream &os, Topology &topo)
  {
    os << "a topology";
    return os;
  }
}

