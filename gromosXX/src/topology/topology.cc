/**
 * @file topology.cc
 * methods definition
 */

namespace simulation
{
  class Multibath;
}

#include <util/stdheader.h>

#include <util/debug.h>
#include <topology/core/core.h>

#include <topology/topology.h>
#include <simulation/multibath.h>

#undef MODULE
#undef SUBMODULE
#define MODULE topology
#define SUBMODULE topology

/**
 * Constructor
 */
topology::Topology::Topology()
  : m_mass(0),
    m_charge(0),
    m_num_solute_chargegroups(0)
{
  m_chargegroup.push_back(0);
  m_molecule.push_back(0);
}

/**
 * set the capacity of solute atoms by resizeing
 * the apropriate arrays.
 */
void topology::Topology::resize(size_t const atoms)
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
void topology::Topology
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
 * add solvent molecules to the simulation (system).
 */
void topology::Topology::solvate(size_t solv, size_t num_molecules)
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
    DEBUG(8, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    m_molecule.push_back(n);

  }
    
}

/**
 * total number of solvent atoms.
 */
size_t topology::Topology::num_solvent_atoms()const
{
  size_t n = 0;
  for(std::vector<size_t>::const_iterator it = m_num_solvent_atoms.begin(),
	to = m_num_solvent_atoms.end();
      it != to; ++it)
    n += *it;
  return n;
}

/**
 * calculate constraint degrees of freedom.
 */
void
topology::Topology::
calculate_constraint_dof(simulation::Multibath &multibath)const
{
  // substract constraints
  {
    std::vector<two_body_term_struct>::const_iterator 
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
  
  DEBUG(7, "and the perturbd distance constraints (DOF calc)");
  
  {
    // substract perturbed constraints
    std::vector<perturbed_two_body_term_struct>
      ::const_iterator 
      c_it = perturbed_solute().distance_constraints().begin(),
      c_to = perturbed_solute().distance_constraints().end();
    
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
  }
  
  DEBUG(10, "end dof calc");

}

void
topology::Topology::update_for_lambda()
{
  for(std::map<size_t, topology::Perturbed_Atom>::const_iterator
        it = perturbed_solute().atoms().begin(),
        to = perturbed_solute().atoms().end();
      it != to; ++it){
    mass()(it->second.sequence_number()) = 
      (1-lambda()) * it->second.A_mass() + 
      lambda() * it->second.B_mass();

    DEBUG(8, "mass A : " << it->second.A_mass() << " B : "
          << it->second.B_mass());
    DEBUG(8, "mass(" << it->second.sequence_number()
          << ") = " << mass()(it->second.sequence_number()));
  }

  // this is nowadays done directly in the shake routines
  // perturbed_solute().set_distance_constraints(lambda());
  
}

/**
 * check state
 */
int 
topology::Topology::check_state()const
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

namespace topology
{
  /**
   * output information about the topology.
   */
  std::ostream & operator<<(std::ostream &os, Topology &topo)
  {
    os << "a topology";
    return os;
  }
}

