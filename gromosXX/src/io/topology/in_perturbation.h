/**
 * @file in_perturbation.h
 * read in a perturbation topology file (03 format)
 */

#ifndef INCLUDED_IN_PERTURBATION_H
#define INCLUDED_IN_PERTURBATION_H

#include <gromosXX/io/instream.h>

namespace io
{
  /**
   * @class In_Perturbation
   * reads in a perturbation topology file (03 version)
   * and parses it into Topology
   * @sa topology::Topology
   */
  class In_Perturbation : public GInStream
  {
  public:
    /**
     * Constructor.
     */
    In_Perturbation(std::istream &is);
    /**
     * parse the topology.
     */
    void read(topology::Topology &topo, simulation::Parameter &param);
    
  };
  
} // io

#endif
