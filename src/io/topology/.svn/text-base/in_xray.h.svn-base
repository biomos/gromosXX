/**
 * @file in_xray.h
 * read in a x-ray restraining file.
 */
/**
 * @page xrayresfile X-ray restraints format
 * @date 12-07-2011
 *
 * A x-ray restraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref xrayresspec
 * - @ref xrayresspec "XRAYRFREESPEC block"
 * - @ref elementmap
 * - @ref elementmap "XRAYSOLVELEMENTSPEC block"
 * - @ref xraybfoccspec
 * - @ref xraybfoccspec "XRAYSOLVBFOCCSPEC block"
 * - @ref xrayrespara
 * - @ref xrayumbrellaweight
 * - @ref xraysymresspec
 * - @ref xraybfopt
 * - @ref xraysfcalc
 * - @ref xrayreplicaexchange
 * - @ref xrayoverallbfactor
 */

#ifndef INCLUDED_IN_XRAYRES_H
#define INCLUDED_IN_XRAYRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Xrayresspec
   * reads in a xray restraining specification file
   */
  class In_Xrayresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Xrayresspec() {}
    /**
     * Constructor.
     */
    In_Xrayresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a xray restraining specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
  
} // io

#endif


