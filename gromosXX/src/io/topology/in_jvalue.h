/**
 * @file in_jvalue.h
 * read in a jvalue restraining specification file
 */
/**
 * @page jvalue J-value restraints specification format
 * @date 28-10-2008
 *
 * A J-value restraining specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref jvalresspec
 */

#ifndef INCLUDED_IN_JVALUE_H
#define INCLUDED_IN_JVALUE_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Jvalue
   * Reads in a @f$^3J@f$-value restraining specification file.
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
     * Read in a @f$^3J@f$-value restraining specification file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
