/**
 * @file dihedral_interaction.h
 * dihedral interaction.
 */

#ifndef INCLUDED_DIHEDRAL_INTERACTION_H
#define INCLUDED_DIHEDRAL_INTERACTION_H

namespace interaction
{
  /**
   * @class Dihedral_Interaction
   * calculates the dihedral interactions.
   */
  template<typename t_interaction_spec>
  class Dihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_Interaction() : Interaction("Dihedral") {}
    /**
     * Destructor.
     */
    virtual ~Dihedral_Interaction() {}

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * the angle type parameters.
     */
    std::vector<dihedral_type_struct> const & parameter()const { return m_parameter; }
    /**
     * the angle type parameters.
     */
    std::vector<dihedral_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<dihedral_type_struct> m_parameter;

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

// template methods
#include "dihedral_interaction.cc"

#endif
