/**
 * @file in_tf_rdc.h
 * read in a tensor-free RDC restraints file.
 */
/**
 * @page tfrdc tensor-free RDC restraints format
 * @date 25-04-2019
 *
 * A tensor-free RDC restraints specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref tfrdcresspec
 */


#ifndef IN_TF_RDC_H
#define	IN_TF_RDC_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Tfrdcresspec
   * reads in a tensor-free RDC restraints file
   */
  class In_Tfrdcresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Tfrdcresspec() {}
    /**
     * Constructor.
     */
    In_Tfrdcresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a tensor-free RDC restraints file.
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

#endif	/* IN_TF_RDC_H */
