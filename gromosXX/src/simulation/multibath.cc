/**
 * @file multibath.cc
 * the multibath parameter class.
 */

#include <stdheader.h>

#include <simulation/multibath.h>
#include <simulation/simulation.h>

#include <configuration/energy.h>
#include <topology/topology.h>

#include <util/error.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE simulation

void simulation::Multibath
::calculate_degrees_of_freedom(topology::Topology &topo,
			       bool rottrans_constraints)
{
  // check whether we have at least one bath
  if (size() == 0){
    io::messages.add("Adding a bath, no temperature coupling",
		     "Multibath::calculate_degrees_of_freedom",
		     io::message::notice);
    add_bath(0.0);
  }
  
  // check whether the last last is really the last_atom
  if (m_bath_index.size() == 0 || 
      (m_bath_index.end()-1)->last_atom != topo.num_atoms()-1){

    io::messages.add("Adding atoms to the last bath!",
		     "Multibath::calculate_degrees_of_freedom",
		     io::message::notice);
    add_bath_index(topo.num_atoms()-1, unsigned(topo.molecules().size())-1, 
		   unsigned(size())-1, unsigned(size())-1);
    
  }

  // loop over the ranges
  std::vector<bath_index_struct>::iterator it = m_bath_index.begin(),
    to = m_bath_index.end();

  DEBUG(8, "number of baths: " << unsigned(size()));
  DEBUG(8, "molecules " << unsigned(topo.molecules().size()));
  
  for(int last=-1; it != to; ++it){
    // get the number of molecules in the range
    int num_mol = 0;
    int mol = 0;

    DEBUG(8, "last atom: " << it->last_atom);
    DEBUG(8, "end of last group: " << last);

    for(topology::Molecule_Iterator m_it = topo.molecule_begin(),
	  m_to = topo.molecule_end();
	m_it != m_to;
	++m_it, ++mol){

      DEBUG(8, "current mol begins: " << (*m_it.begin()));
      
      if ((*(m_it.begin())) > it->last_atom){
	break;
      }
      
      if (int(*(m_it.begin())) > last)
	++num_mol;
    }

    // add the last molecule
    DEBUG(8, "last molecule is " << mol - 1);
    it->last_molecule = mol - 1;

    DEBUG(8, "num molecules is " << num_mol);
    // add the molecular translational dof
    (*this)[it->com_bath].dof += num_mol * 3;
    (*this)[it->com_bath].com_dof += num_mol * 3;

    // and the internal and molecular rotational dof
    (*this)[it->ir_bath].dof += (it->last_atom - last) * 3 - num_mol * 3;
    (*this)[it->ir_bath].ir_dof += (it->last_atom - last) * 3 - num_mol * 3;

    last = it->last_atom;
  }

  // substract constraints
  topo.calculate_constraint_dof(*this, rottrans_constraints);
  
}

int simulation::Multibath::check_state(unsigned int num_atoms)const
{
  int result = 0;
  unsigned int last_atom = 0;
  std::vector<bath_index_struct>::const_iterator it = m_bath_index.begin(),
    to = m_bath_index.end();
  for( ; it!=to; ++it){
    if (it->last_atom < last_atom){
      io::messages.add("Multibath not sorted", "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if (it->last_atom > num_atoms){
      io::messages.add("Multibath last atom index too large",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if (it->com_bath >= size()){
      io::messages.add("Multibath: com bath index out of range",
		       "Multibath::check_state",
		       io::message::error);
      return E_INPUT_ERROR;
    }
    if (it->ir_bath >= size()){
      io::messages.add("Multibath: ir bath index out of range",
		       "Multibath::check_state",
		       io::message::error);
      return E_INPUT_ERROR;
    }    
    if ((*this)[it->com_bath].dof == 0)
      io::messages.add("Multibath: bath with 0 degrees of freedom?",
		       "Multibath::check_state",
		       io::message::warning);
    if ((*this)[it->ir_bath].dof == 0)
      io::messages.add("Multibath: bath with 0 degrees of freedom?",
		       "Multibath::check_state",
		       io::message::warning);
    if ((*this)[it->ir_bath].solute_constr_dof < 0 
	|| (*this)[it->ir_bath].solvent_constr_dof < 0){
      io::messages.add("Multibath: constrained degrees of freedom negative",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if ((*this)[it->com_bath].solute_constr_dof < 0 
	|| (*this)[it->com_bath].solvent_constr_dof < 0){
      io::messages.add("Multibath: constrained degrees of freedom negative",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if ((*this)[it->ir_bath].tau < 0 && (*this)[it->ir_bath].tau != -1){
      io::messages.add("Multibath: tau < 0 && tau != -1",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if ((*this)[it->com_bath].tau < 0 && (*this)[it->com_bath].tau != -1){
      io::messages.add("Multibath: tau < 0 && tau != -1",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    if ((*this)[it->ir_bath].temperature < 0 ||
	(*this)[it->com_bath].temperature < 0){
      io::messages.add("Multibath: temperature < 0",
		       "Multibath::check_state",
		       io::message::error);
      ++result;
    }
    
  }
  
  return result;

}

void simulation::Multibath::calc_totals(configuration::Energy const &energy,
					double & ekin, double & ekin_mol, 
					double & ekin_ir,
					double & temp, double & temp_mol, 
					double & temp_ir,
					double & scale)const
{
  ekin = 0.0;
  ekin_mol = 0.0;
  ekin_ir = 0.0;
  scale = 0.0;
  
  double sum_dof = 0.0;
  double sum_dof_com = 0.0;
  double sum_dof_ir = 0.0;
  double tau_dof = 0.0;
  
  std::vector<bath_struct>::const_iterator
    it = begin(),
    to = end();
  
  for(unsigned int i=0; it != to; ++it, ++i){
      
    const double e_kin = energy.kinetic_energy[i];
    const double e_kin_com = energy.com_kinetic_energy[i];
    const double e_kin_ir = energy.ir_kinetic_energy[i];

    if (it->tau != -1){
      tau_dof += it->dof;
      scale += it->scale * it->dof;
    }

    sum_dof += it->dof;
    sum_dof_ir += it->ir_dof;
    sum_dof_com += it->com_dof;

    ekin += e_kin;
    ekin_mol += e_kin_com;
    ekin_ir += e_kin_ir;
  }

  if (sum_dof)
    temp = 2 * ekin / (math::k_Boltzmann * sum_dof);
  else temp = 0.0;
  
  if (sum_dof_com)
    temp_mol = 2 * ekin_mol / (math::k_Boltzmann * sum_dof_com);
  else temp_mol = 0.0;
  
  if (sum_dof_ir)
    temp_ir = 2 * ekin_ir / (math::k_Boltzmann * sum_dof_ir);
  else temp_ir = 0.0;
  
  if (tau_dof)
    scale /= tau_dof;
  else scale = 0.0;

}
