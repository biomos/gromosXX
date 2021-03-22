/**
 * @file   in_bsleus.h
 * Read in the file defining the B&S-LEUS Umbrella
 */
/**
 * @page bsleus B&S-LEUS format
 * 
 * A B&S-LEUS topology file may have the following blocks:
 *  - @ref title
 *  - @ref bsleussub
 *  - @ref bsleuscoord
 *  - @ref bsleussph
 *  - @ref bsleusstk
 * 
 * @sa @ref bsleusparam
 * @sa @ref bsleusmem
 * @sa @ref bsleuspos
 * 
 * 
 */

#ifndef IN_BSLEUS_H
#define	IN_BSLEUS_H

#include "../instream.h"

namespace io {
  /**
   * @class In_BSLEUS reads in the definition of the B&S-LEUS Umbrella 
   * definition.
   */
  class In_BSLEUS : public GInStream {
  public:
    /**
     * Constructor
     * @param is the file to read in.
     */
    In_BSLEUS (std::istream &is) : GInStream(is) {readStream();}
    /**
     * Read in the definition file of the B&S-LEUS algorithm
     * @param topo the topology
     * @param sim the simulation parameter
     * @param os  
     */
    void read(topology::Topology &topo,
          configuration::Configuration &conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  private:
    /**
     * Find out, wheter there was an error in getting the tokens
     * @param pos
     */
    void findError(size_t pos);
    /**
     * Parse the Atomspecifier in coordStr
     * @param      conf
     * @param[in]  coordStr The atom specifier (or a simple number)
     * @param[in]  refFiles The reference Files
     * @param[out] coords   The coordinates
     */
    void parseSpecifier(topology::Topology &topo,
                        simulation::Simulation &sim,
                        configuration::Configuration &conf,
                        std::string &coordStr,
                        std::map<unsigned int, std::string> &refFiles,
                        std::vector<unsigned int> &cartAtoms,
                        std::vector<double> &coords,
                        std::ostream & os);
    /**
     * Put the positions into a box centered at zero
     * @param conf
     * @param pos
     */
    void put_into_box(configuration::Configuration &conf, math::VArray& pos);

    template<math::boundary_enum B>
    void _put_into_box(configuration::Configuration& conf, math::VArray& pos);
  };
}
#endif	/* IN_BSLEUS_H */

