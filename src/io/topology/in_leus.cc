/**
 * @file in_leus.cc
 * implements methods of In_Localelevspec
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../simulation/parameter.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/topology/in_leus.h"
#include "../../io/configuration/in_configuration.h"
#include "../../util/le_coordinate.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section localelevspec LOCALELEVSPEC block
 * The LOCALELEVSPEC specifies the dihedral angles which we attach to
 * umbrella potentials.
 *
 * The block is read from the local elevation specification file
 * (\@led).
 *
 * @verbatim
LOCALELEVSPEC
#  NLEPID  IPLE  JPLE  KPLE  LPLE
        1     1     2     3     4
        1     2     3     4     5
        2     1     2     3     4
END
@endverbatim
 */
void
io::In_Localelevspec::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a local elevation definition file");

  std::vector<std::string> buffer;

  { // LOCALELEVSPEC

    buffer = m_block["LOCALELEVSPEC"];
    DEBUG(10, "LOCALELEVSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no LOCALELEVSPEC block in local elevation definition file",
		       "In_Localelevspec", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in LOCALELEVSPEC data");

    int id;
    unsigned int n, i, j, k, l;
    for(n=0; it != to; ++n, ++it){
      DEBUG(11, "\tnr " << n);

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> id >> i >> j >> k >> l;
      if (_lineStream.fail()){
        io::messages.add("bad line in LOCALELEVSPEC block",
                "In_Localelevspec", io::message::error);
        return;
      }

      // convert to gromos
      --i; --j; --k; --l;

      if (i > topo.num_atoms() ||
          j > topo.num_atoms() ||
          k > topo.num_atoms() ||
          l > topo.num_atoms()) {
        std::ostringstream msg;
        msg << "Dihedral (" << i+1 << "-" << j+1 << "-" << k+1 << "-" << l+1
                << ") atom indices out of range.";
        io::messages.add(msg.str(), "In_Localelevspec", io::message::error);
        return;
      }

      util::LE_Dihedral_Coordinate * dih = new util::LE_Dihedral_Coordinate(id, i, j, k, l);
      topo.le_coordinates().push_back(dih);
    }
  } // LOCALELEVSPEC
}

/**
 * @section leusbias LEUSBIAS block
 * The LOCALELEVSPEC specifies the dihedral angles which we attach to
 * umbrella potentials.
 *
 * The block is read from the local elevation umbrella sampling database file
 * (\@lud) or the configuration. The reading from the input file is not implemented
 * in md++.
 *
 * @verbatim
LEUSBIAS
# NUMB:    number of umbrellas
# NLEPID:  ID of the umbrella
# NDIM:    number of dimensions of the umbrella
# CLES:    force constant of the umbrella
# VARTYPE: type of the collective variable
# NTLEFU:  functional form of the potential in this dimension
# WLES:    width of the potential in grid units
# RLES:    cutoff in grid units
# NGRID:   number of grid points
# GRIDMIN: minumum value for the grid (origin)
# GRIDMAX: maximum value for the grid.
#    if GRIDMIN == GRIDMAX: periodic coordinate such as angle
# NCONLE:  number of (different) configurations visited
# NVISLE:  number of times the confiuguration was visited
# ICONF(1..NDIM): the coordinates of the grid point (1..GRIDMAX)
#
# NUMB
2
# NLEPID NDIM CLES
         3         1    0.005000000
# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)
         1         0    1.000000000    1.000000000        32    0.000000000    0.000000000
# NCONLE
         3
# NVISLE ICONF(1..NDIM)
         2         3
         5         4
         4         6
# NLEPID NDIM CLES
         7         2    0.004000000
# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)
         1         1    1.000000000    1.000000000        16   15.000000000   15.000000000
         1         0    1.000000000    1.000000000        10    0.000000000    0.000000000
# NCONLE
         1
# NVISLE ICONF(1..NDIM)
         1         3         7
END
@endverbatim
 */
void
io::In_LEUSBias::read(topology::Topology& topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os) {

  DEBUG(7, "reading in a local elevation definition file");

  std::vector<std::string> buffer;

  if (sim.param().localelev.read) { // LOCALELEVSPEC
    buffer = m_block["LEUSBIAS"];
    DEBUG(10, "LOCALELEVSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no LEUSBIAS block in LEUS database file",
		       "In_Localelevspec", io::message::error);
      return;
    }

    // remove the title from the buffer. the title to the buffer
    std::vector<std::string> new_buffer;
    new_buffer.insert(new_buffer.begin(), buffer.begin()+1, buffer.end());

    io::In_Configuration::_read_leusbias(conf.special().umbrellas, new_buffer, false);

  } // LEUSBIAS
}


