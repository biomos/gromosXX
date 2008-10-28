/**
 * @file in_jvalue.h
 * read in a jvalue restraining specification file.
 */

#ifndef INCLUDED_IN_JVALUE_H
#define INCLUDED_IN_JVALUE_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Jvalue
   * reads in a jvalue restraining specification file.
   */
  class In_Jvalue : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Jvalue() {}
    /**
     * Constructor.
     */
    In_Jvalue(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraining file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
