/**
 * @file in_posres.h
 * read in a position restraining file.
 */

#ifndef INCLUDED_IN_POSRES_H
#define INCLUDED_IN_POSRES_H

namespace io {

  /**
   * @class In_Topology
   * reads in a position restraining file
   */
  class In_Posres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Posres() {}
    /**
     * Constructor.
     */
    In_Posres(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraining file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim);

  };

} // io

#endif