/**
 * @file multibath.tcc
 * methods for the multibath parameter class.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE simulation

#include "../../debug.h"

inline
simulation::Multibath::Multibath()
{
}

inline void simulation::Multibath
::add_bath(simulation::bath_struct s)
{
  push_back(s);
}

inline void simulation::Multibath
::add_bath(double temperature, double tau,
	   double dof, double com_dof,
	   double ir_dof, 
	   double solute_constr_dof, double solvent_constr_dof)
{
  bath_struct s;
  s.temperature = temperature;
  s.tau = tau;
  s.dof = dof;
  s.ir_dof = ir_dof;
  s.com_dof = com_dof;
  s.solute_constr_dof = solute_constr_dof;
  s.solvent_constr_dof = solvent_constr_dof;
  
  add_bath(s);
}

inline void simulation::Multibath
::add_bath_index(bath_index_struct s)
{
  m_bath_index.push_back(s);
}

inline void simulation::Multibath
::add_bath_index(size_t const last, size_t const com_bath, size_t const ir_bath)
{
  bath_index_struct s;
  s.last_atom = last;
  s.com_bath = com_bath;
  s.ir_bath = ir_bath;

  add_bath_index(s);
}

inline simulation::bath_struct & 
simulation::Multibath::bath(size_t i)
{
  assert(i < size());
  return (*this)[i];  
}

inline simulation::bath_struct const & 
simulation::Multibath::bath(size_t i)const
{
  assert(i < size());
  return (*this)[i];  
}

inline std::vector<simulation::bath_index_struct> & 
simulation::Multibath::bath_index()
{
  return m_bath_index;
}

inline std::vector<simulation::bath_index_struct> const & 
simulation::Multibath::bath_index()const
{
  return m_bath_index;
}

inline void simulation::Multibath::in_bath(size_t const atom,
					   size_t &com, size_t &ir)const
{
  std::vector<bath_index_struct>::const_iterator it = m_bath_index.begin(),
    to = m_bath_index.end();
  
  for(; it != to; ++it){
    if (it->last_atom >= atom){
      com = it->com_bath;
      ir = it->ir_bath;
      return;
    }
  }

  // if no bath read in, the 0 bath is not yet added (change that!)
  assert(false);

  com = 0;
  ir = 0;  
}

template<typename t_topology>
inline void simulation::Multibath
::calculate_degrees_of_freedom(t_topology &topo)
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
    add_bath_index(topo.num_atoms()-1, size()-1, size()-1);
    
  }

  // loop over the ranges
  std::vector<bath_index_struct>::iterator it = m_bath_index.begin(),
    to = m_bath_index.end();

  for(int last=-1; it != to; ++it){
    // get the number of molecules in the range
    int num_mol = 0;
    int mol = 0;
    for(Molecule_Iterator m_it = topo.molecule_begin(),
	  m_to = topo.molecule_end();
	m_it != m_to;
	++m_it, ++mol){
      if (*(m_it.begin()) > it->last_atom){
	break;
      }
      
      if (int(*(m_it.begin())) > last)
	++num_mol;
    }
    // add the last molecule
    it->last_molecule = mol - 1;
    // add the molecular translational dof
    (*this)[it->com_bath].dof += num_mol * 3;
    (*this)[it->com_bath].com_dof += num_mol * 3;

    // and the internal and molecular rotational dof
    (*this)[it->ir_bath].dof += (it->last_atom - last) * 3 - num_mol * 3;
    (*this)[it->ir_bath].ir_dof += (it->last_atom - last) * 3 - num_mol * 3;

    last = it->last_atom;
  }

  // substract constraints
  topo.calculate_constraint_dof(*this);
  
}

inline int simulation::Multibath::check_state(size_t const num_atoms)const
{
  int result = 0;
  size_t last_atom = 0;
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
      throw std::string("com bath index out of range");
    }
    if (it->ir_bath >= size()){
      io::messages.add("Multibath: ir bath index out of range",
		       "Multibath::check_state",
		       io::message::error);
      throw std::string("ir bath index out of range");
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
