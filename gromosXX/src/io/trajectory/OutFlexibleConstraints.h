/**
 * @file OutFlexibleConstraints.h
 * write out flexible constraint information.
 */

#ifndef INCLUDED_OUTFLEXIBLECONSTRAINTS_H
#define INCLUDED_OUTFLEXIBLECONSTRAINTS_H

namespace io
{
  /**
   * @class OutFlexibleConstraints
   * print out flexible constraints information
   */
  class OutFlexibleConstraints
  {
  public:
    /**
     * Constructor.
     */
    OutFlexibleConstraints(std::ostream & os);
    

    /**
     * write out a flexible constraints file from a topology.
     */
    void write_FLEXCON(std::vector<double> &vel,
		       simulation::Topology &topo);

    /**
     * write out a flexible constraints file from a perturbation topology.
     */
    void write_FLEXCON(std::vector<double> &vel,
		       simulation::Perturbation_Topology &topo);

    /**
     * write a title block.
     */
    void write_title(std::string title);

  protected:
    std::ostream & m_os;
  };
} // io

#include "OutFlexibleConstraints.tcc"

#endif


  
