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
     * Read in a G96 trajectory into the system.
     */
    void read_FLEXCON(std::vector<double> &vel,
		      std::vector<simulation::compound::
		      distance_constraint_struct> &constr);
  private:
    
  };
  
  
} // io

// template and inline methods
#include "InFlexibleConstraints.tcc"

#endif
