/**
 * @file harmonic_bond_interaction.h
 * harmonic bond interaction.
 */

#ifndef INCLUDED_HARMONIC_BOND_INTERACTION_H
#define INCLUDED_HARMONIC_BOND_INTERACTION_H

namespace interaction
{
  /**
   * @class harmonic_bond_interaction
   * calculates the bond interactions (harmonic).
   */
  template<typename t_simulation>
  class harmonic_bond_interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Destructor.
     */
    virtual ~harmonic_bond_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);
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
#include "harmonic_bond_interaction.tcc"

#endif
