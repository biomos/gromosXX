/**
 * @file dihedral_interaction.h
 * dihedral interaction.
 */

#ifndef INCLUDED_DIHEDRAL_INTERACTION_H
#define INCLUDED_DIHEDRAL_INTERACTION_H

namespace interaction
{
  /**
   * @class Dihedral_interaction
   * calculates the dihedral interactions.
   */
  template<typename t_simulation>
  class Dihedral_interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_interaction();
    /**
     * Destructor.
     */
    virtual ~Dihedral_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);
    /**
     * add improper dihedral type.
     */
    void add(dihedral_type_struct s);
    /**
     * add improper type.
     */
    void add(double K, double pd, int m);
    /**
     * the angle type parameters.
     */
    std::vector<dihedral_type_struct> const & parameter()const;
    
  protected:
    std::vector<dihedral_type_struct> m_dihedral_parameter;
    
  };
  
} // interaction

// template methods
#include "dihedral_interaction.tcc"

#endif
