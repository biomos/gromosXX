/**
 * @file improper_dihedral_interaction.h
 * improper dihedral interaction.
 */

#ifndef INCLUDED_IMPROPER_DIHEDRAL_INTERACTION_H
#define INCLUDED_IMPROPER_DIHEDRAL_INTERACTION_H

namespace interaction
{
  /**
   * @class Improper_dihedral_interaction
   * calculates the improper dihedral interactions.
   */
  template<typename t_simulation>
  class Improper_dihedral_interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Improper_dihedral_interaction();
    /**
     * Destructor.
     */
    virtual ~Improper_dihedral_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);
    /**
     * add improper dihedral type.
     */
    void add(improper_dihedral_type_struct s);
    /**
     * add improper type.
     */
    void add(double K, double q);
    /**
     * the angle type parameters.
     */
    std::vector<improper_dihedral_type_struct> const & parameter()const;
    
  protected:
    std::vector<improper_dihedral_type_struct> m_improper_dihedral_parameter;
    
  };
  
} // interaction

// template methods
#include "improper_dihedral_interaction.tcc"

#endif
