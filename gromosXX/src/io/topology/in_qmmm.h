/**
 * @file in_qmmm.h
 * read in a QM/MM specification file
 */
/**
 * @page qmmm QM/MM specification format
 * @date 22-01-2020
 *
 * A QM/MM specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref qmzone
 * - @ref qmunit
 */

#ifndef INCLUDED_IN_QMMM_H
#define INCLUDED_IN_QMMM_H

#include "../instream.h"

namespace io {

  /**
   * @class In_QMMM
   * reads in a QM/MM specification file
   */
  class In_QMMM : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_QMMM() {}
    /**
     * Constructor.
     */
    In_QMMM(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a QM/MM specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Read in QM/MM units conversion factors
     */
    void read_units(const simulation::Simulation & sim
          , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read the map of atomic numbers to element names
     */
    void read_elements(const topology::Topology& topo
    , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param);
    /**
     * Read the list of QM atoms
     */
    void read_zone(topology::Topology& topo
                    , simulation::Simulation& sim
                    , const std::string& blockname);
  };
} // io

#endif
