/**
 * @file quartic_bond_interaction.h
 * quartic bond interaction.
 */

#ifndef INCLUDED_QUARTIC_BOND_INTERACTION_H
#define INCLUDED_QUARTIC_BOND_INTERACTION_H

namespace interaction
{
  /**
   * @class Quartic_bond_interaction
   * calculates the bond interactions (quartic).
   */
  template<typename t_simulation, typename t_interaction_spec>
  class Quartic_bond_interaction : 
    public Interaction<t_simulation, t_interaction_spec>
  {
  public:
    /**
     * Constructor.
     */
    Quartic_bond_interaction();
    /**
     * Destructor.
     */
    virtual ~Quartic_bond_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);
    /**
     * add bond type.
     */
    void add(bond_type_struct s);
    /**
     * add bond type.
     */
    void add(double K, double r0);
    /**
     * the bond type parameters.
     */
    std::vector<bond_type_struct> const & parameter()const;
    
  protected:
    std::vector<bond_type_struct> m_bond_parameter;
    
  };
  
} // interaction

// template methods
#include "quartic_bond_interaction.tcc"

#endif
