/**
 * @file:   in_reference.h
 * Interface for reading in a refrence
 */

#ifndef IN_REFERENCE_H
#define	IN_REFERENCE_H

#include <vector>
#include "../instream.h"

namespace math {
  class VArray;
}

namespace io {
  /**
   * @class In_Reference
   * Read in a reference structure and to which atoms, one should fit.
   */
  class In_Reference : public GInStream {
  public:
    /**
     * Default Constructor
     */
    In_Reference() {}
    /**
     * Constructor
     * @param is    file stream to read in
     */
    In_Reference(std::istream& is) : GInStream(is) {
      readStream();
    }
    /**
     * Read in the file with a reference positions.
     */
    void read(topology::Topology &topo,
                simulation::Simulation &sim, 
                configuration::Configuration& conf,
                std::ostream &os = std::cout);
    math::VArray refpos(){ return m_refpos;}
  private:
    math::VArray m_refpos;
  };
}


#endif	/* IN_REFERENCE_H */

