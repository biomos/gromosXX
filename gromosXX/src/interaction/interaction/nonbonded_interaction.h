/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_simulation, typename t_pairlist>
  class Nonbonded_Interaction : public Interaction<t_simulation>
  {
  public:    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

    /**
     * add the lj parameters for atom type i and j.
     */
    void add_lj_parameter(size_t iac_i, size_t iac_j,
			  lj_parameter_struct lj);
    
    /**
     * get the lj parameters for atom type i and j.
     */
    lj_parameter_struct const & lj_parameter(size_t iac_i, size_t iac_j);
    
    /**
     * resize the lj_parameter matrix.
     */
    void resize(size_t i);

    /**
     * pairlist accessor
     */
    t_pairlist & pairlist();
    
  protected:
    /**
     * helper class to build the pairlist.
     * @TODO parametrize that one?
     * or assume an iterator and take a reference (polymorphism)
     */
    t_pairlist m_pairlist;

    /**
     * the lj parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

    /**
     * the long-range force.
     */
    math::VArray m_longrange_force;
    
    /**
     * calculate the interactions.
     */
    void do_interactions(t_simulation &sim,
			 typename t_pairlist::iterator it, 
			 typename t_pairlist::iterator to,
			 math::VArray &force);

    /**
     * calculate the 1,4-interactions.
     */
    void do_14_interactions(typename t_simulation::topology_type &topo,
			    typename t_simulation::system_type &sys);

  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
