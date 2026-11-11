/**
 * @file in_colvarres.h
 * read in a collective variable restraints file.
 */
/**
 * @page colvarres colvar restraints format
 * @date 13-11-2025
 *
 * A colvar restraints specification file may contain the following blocks:
 * - @ref title
 * - @ref coordnumresspec
 */

#ifndef INCLUDED_IN_COLVARRES_H
#define INCLUDED_IN_COLVARRES_H

#include "../instream.h"
#include <set>
#include <string>

namespace io {

  /**
   * @class In_Colvarres
   * reads in a collective variable restraints file
   */
  class In_Colvarres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Colvarres() {}
    /**
     * Constructor.
     */
    In_Colvarres(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a colvar restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;

  private:
    // Per-reader bookkeeping: no state is shared with other input readers.
    std::set<std::string> m_colvar_blocks_read;
    void read_restraints(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_pertrestraints(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_distanceresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_angresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_dihedralresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_coordnumresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_pertdistanceresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_pertangresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_pertdihedralresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
    void read_pertcoordnumresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os);
  };

} // io

#endif