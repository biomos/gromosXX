/**
 * @file in_reference.cc
 * Implements methods from the class In_Reference
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/configuration/in_configuration.h"
#include "../../io/configuration/out_configuration.h"

#include "../../math/gmath.h"
#include "../../math/transformation.h"
#include "in_reference.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology


/**
 * @section refconf REFCONF block
 * 
 * Defines the reference configuration. It is read from a file, whose filename
 * is defined in the bsleus topology file in the @ref bsleuscoord.
 * 
 * @verbatim
REFCONF
    1 ALA   CA         1    0.099562952    0.066728937    0.299959523
    1 ALA   C          2    0.033282983    0.063516406    0.162089587
    1 ALA   O          3    0.033835294    0.163010157    0.089764673
    1 ALA   N          4   -0.019654805   -0.053770863    0.128469837
    1 ALA   H          5   -0.010672571   -0.132310509    0.189720041
    1 ALA   CA         6   -0.090279524   -0.084088083    0.003152546
    1 ALA   CB         7   -0.220622683   -0.004933744   -0.009448523
    1 ALA   C          8   -0.007031111   -0.075815957   -0.124950045
    1 ALA   O          9    0.033835294   -0.179846848   -0.176319824
    1 ALA   N         10    0.026861616    0.046083612   -0.165973935
    1 ALA   H         11    0.021319603    0.124697956   -0.104405500
    1 ALA   CA        12    0.099562952    0.066728937   -0.292058378
    .....
END
@endverbatim
 * 
 * @sa @ref bsleuscoord
 */
void
io::In_Reference::read(topology::Topology& topo, 
        simulation::Simulation& sim, 
        configuration::Configuration& conf,
        std::ostream& os)
{
  DEBUG(7, "reading in a reference file");

  if (!quiet)
    os << "Reading in a Reference configuration\n";
  
  std::vector<std::string> buffer;
  
  // read the reference positions from this file
  { // REFPOSITION
    buffer = m_block["POSITION"];

    if (!buffer.size()) {
      io::messages.add("no POSITION block in reference position file!", 
              "In_reference", io::message::error);
      return;
    } else {
      // remove title line.
      DEBUG(8, "Read in the reference positions.");
      m_refpos.resize(topo.num_atoms());
      buffer.erase(buffer.begin(), buffer.begin()+1);
      io::In_Configuration::_read_position(m_refpos, buffer, topo.num_atoms(), 
                                           std::string("POSITION"));
      if (conf.boundary_type == math::truncoct) {
        math::truncoct_triclinic(m_refpos, true);
      }
    }
  } // REFPOSITION
}
