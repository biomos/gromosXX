/**
 * @file multibath.h
 * the multibath parameter class.
 */

#ifndef INCLUDED_MULTIBATH_H
#define INCLUDED_MULTIBATH_H

namespace simulation
{
  /**
   * @struct bath_struct
   * holds the bath / degree of freedom information
   */
  struct bath_struct
  {
    double temperature;
    double tau;
    double dof;
    double ir_dof;
    double com_dof;
    double solute_constr_dof;
    double solvent_constr_dof;
    double scale;
  };

  /**
   * @struct bath_index_struct
   * holds bath index for a range of atoms.
   */
  struct bath_index_struct
  {
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
    Multibath();

    /**
     * add a bath.
     */
    void add_bath(bath_struct s);
    /**
     * add a bath.
     */
    void add_bath(double temperature,
		  double tau = -1, double dof = 0, 
		  double com_dof = 0, double ir_dof = 0,
		  double solute_constr_dof = 0, double solvent_constr_dof = 0);

    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(bath_index_struct s);

    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(size_t const last, size_t const com_bath, size_t const ir_bath);
    
    /**
     * get bath i.
     */
    bath_struct & bath(size_t i);
    /**
     * get const bath i.
     */
    bath_struct const & bath(size_t i)const;
    /**
     * bath indices accessor.
     */
    std::vector<bath_index_struct> & bath_index();
    /**
     * const bath indices accessor.
     */
    std::vector<bath_index_struct> const & bath_index()const;
    /**
     * get the bath number of particle number i.
     */
    void in_bath(size_t const atom,
		 size_t &com, size_t &ir)const;

    /**
     * calculate degrees of freedom.
     */
    template<typename t_topology>
    void calculate_degrees_of_freedom(t_topology &topo);
    
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

#include "multibath.tcc"

#endif
