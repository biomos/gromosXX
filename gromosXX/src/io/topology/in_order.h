/**
 * @file in_order.h
 * read in a order parameter restraints file.
 */
/**
 * @page orderparam order parameter restraints format
 * @date 12-04-2011
 *
 * A distance restraints specifcation file may contain the following blocks:
 * - @ref title
 * - @ref orderparamresspec
 */


#ifndef IN_ORDER_H
#define	IN_ORDER_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Orderparamresspec
   * reads in a order parameter restraints file
   */
  class In_Orderparamresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Orderparamresspec() {}
    /**
     * Constructor.
     */
    In_Orderparamresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a order parameter restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif	/* IN_ORDER_H */

