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
    size_t last_atom;
    double temperature;
    double tau;
    double dof;
    double solute_constr_dof;
    double solvent_constr_dof;
    double kinetic_energy;
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
    void add_bath(int last_atom, double temperature,
		  double tau = -1, double dof = 0,
		  double solute_constr_dof = 0, double solvent_constr_dof = 0);

    /**
     * get bath i.
     */
    bath_struct & bath(size_t i);
    /**
     * get const bath i.
     */
    bath_struct const & bath(size_t i)const;
    /**
     * get the bath number of particle number i.
     */
    int in_bath(size_t const i)const;

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
    
  };
  
} // simulation

#include "multibath.tcc"

#endif
