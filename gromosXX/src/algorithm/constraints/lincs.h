/**
 * @file lincs.h
 * the LINCS algorithm.
 */

#ifndef INCLUDED_LINCS_H
#define INCLUDED_LINCS_H

namespace algorithm
{
  /**
   * @class Lincs
   * implements the lincs algorithm.
   */
  template<math::virial_enum do_virial>
  class Lincs : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Lincs();
    
    /**
     * Destructor.
     */
    virtual ~Lincs();
        
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }

    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }

    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim,
	     bool quiet = false);

    virtual void print_timing(std::ostream & os);

  protected:

    std::vector<interaction::bond_type_struct> m_parameter;

    double m_solvent_timing;

  };
  
} //algorithm

// template methods
#include "lincs.tcc"

#endif
