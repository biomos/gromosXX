/**
 * @file multibath.h
 * the multibath parameter class.
 */

#ifndef INCLUDED_MULTIBATH_H
#define INCLUDED_MULTIBATH_H

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE simulation

#include <util/debug.h>

namespace simulation
{
  /**
   * @struct bath_struct
   * holds the bath / degree of freedom information
   */
  struct bath_struct
  {
    bath_struct(double t, double tau, double dof, double ir_dof, 
		double com_dof, double solu_constr_dof, 
		double solv_constr_dof, double scale=0, double ekin=0)
      : temperature(t),tau(tau), dof(dof), ir_dof(ir_dof),com_dof(com_dof),
	solute_constr_dof(solu_constr_dof), 
	solvent_constr_dof(solv_constr_dof), scale(scale), ekin(ekin)
    {}
    
    double temperature;
    double tau;
    double dof;
    double ir_dof;
    double com_dof;
    double solute_constr_dof;
    double solvent_constr_dof;
    double scale;
    double ekin;
  };

  /**
   * @struct bath_index_struct
   * holds bath index for a range of atoms.
   */
  struct bath_index_struct
  {
    bath_index_struct(size_t last_atom, size_t last_molecule, 
		      size_t com_bath, size_t ir_bath)
      : last_atom(last_atom), last_molecule(last_molecule), 
	com_bath(com_bath), ir_bath(ir_bath){}
    
    size_t last_atom;
    size_t last_molecule;
    size_t com_bath;
    size_t ir_bath;
  };
  
  /**
   * @class Multibath
   * holds multibath and degree of freedom information.
   */
  class Multibath : public std::vector<bath_struct>
  {
  public:
    
    /**
     * Constructor.
     */
    Multibath(){};

    /**
     * add a bath.
     */
    void add_bath(bath_struct s) { push_back(s);}

    /**
     * add a bath.
     */
    void add_bath(double temperature,
		  double tau = -1, double dof = 0, 
		  double com_dof = 0, double ir_dof = 0,
		  double solute_constr_dof = 0, double solvent_constr_dof = 0)
    {
      push_back(bath_struct(temperature, tau, dof, com_dof, ir_dof, 
			    solute_constr_dof, solvent_constr_dof));
    }
      
    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(bath_index_struct s){ m_bath_index.push_back(s);}

    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(size_t const last, size_t const last_m, 
			size_t const com_bath, size_t const ir_bath){
      m_bath_index.push_back(bath_index_struct(last, last_m, 
					       com_bath, ir_bath));
    }
    
    
    /**
     * get bath i.
     */
    bath_struct & bath(size_t i) {
      assert(i < size());
      return (*this)[i];  
    }
    /**
     * get const bath i.
     */
    bath_struct const & bath(size_t i)const{
      assert(i < size());
      return (*this)[i];  
    }
    /**
     * bath indices accessor.
     */
    std::vector<bath_index_struct> & bath_index(){return m_bath_index;}
    
    /**
     * const bath indices accessor.
     */
    std::vector<bath_index_struct> const & bath_index()const{
      return m_bath_index;}
    /**
     * get the bath number of particle number i.
     */
    void in_bath(size_t const atom,
		 size_t &com, size_t &ir)const{
      std::vector<bath_index_struct>::const_iterator 
	it = m_bath_index.begin(),
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
    
    
    /**
     * calculate degrees of freedom.
     */
    void calculate_degrees_of_freedom(topology::Topology &topo);
    
    /**
     * check the state.
     */
    int check_state(size_t const num_atoms)const;

  private:
    /**
     * the bath index for a range of atoms.
     */
    std::vector<bath_index_struct> m_bath_index;
    
  };
  
} // simulation


inline void simulation::Multibath
::calculate_degrees_of_freedom(topology::Topology &topo)
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
    add_bath_index(topo.num_atoms()-1, topo.molecules().size()-1, 
		   size()-1, size()-1);
    
  }

  // loop over the ranges
  std::vector<bath_index_struct>::iterator it = m_bath_index.begin(),
    to = m_bath_index.end();

  DEBUG(8, "number of baths: " << size());
  DEBUG(8, "molecules " << topo.molecules().size());
  
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

#endif
