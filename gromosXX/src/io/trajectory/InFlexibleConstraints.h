/**
 * @file InFlexibleConstraints.h
 * read in flexible constraint information.
 */

#ifndef INCLUDED_INFLEXIBLECONSTRAINTS_H
#define INCLUDED_INFLEXIBLECONSTRAINTS_H

namespace io {

  /**
   * @class InFlexibleConstraints
   * reads in flexible constraints information
   * (distance and velocity)
   * into algorithm::constraints::Flexible_Constraints
   */
  class InFlexibleConstraints : public GInStream {

  public:
    /**
     * Constructor.
     */
    InFlexibleConstraints(std::istream& is);
    /**
     * Read in a flexible constraints file into a topology.
     */
    void read_FLEXCON(std::vector<double> &vel,
		      simulation::Topology &topo);

    /**
     * Read in a flexible constraints file into a perturbation topology.
     */
    void read_FLEXCON(std::vector<double> &vel,
		      simulation::Perturbation_Topology &topo);
    

  private:
    
  };
  
  
} // io

// template and inline methods
#include "InFlexibleConstraints.tcc"

#endif
