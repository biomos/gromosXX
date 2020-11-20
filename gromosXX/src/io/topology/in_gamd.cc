/**
 * @file in_gamd.cc
 * implements methods of In_GAMD
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../simulation/parameter.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/topology/in_gamd.h"
#include "../../io/configuration/in_configuration.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section gamdatoms GAMDATOMS block
 * The GAMDATOMS specifies the atoms that should be accelerated and to which acceleration group they should be asigned.
 *
 * The block is read from the gamd specification file
 * (\@gamd).
 *
 * @verbatim
GAMDATOMS
#  INATOM   FINATOM   AGROUP
   1        230       1
   300      400       1
   500      800       2   

END
@endverbatim
 */
void
io::In_GAMD::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a GAMD definition file");

  std::vector<std::string> buffer;

  { // GAMDATOMS

    buffer = m_block["GAMDATOMS"];
    DEBUG(10, "GAMDATOMS block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no GAMDATOMS block in GAMD definition file",
		       "In_GAMD", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in GAMDATOMS data");

    unsigned int n, iatom, fatom, agroup;
    for(n=0; it != to; ++n, ++it){
        DEBUG(11, "\tnr " << n);

        _lineStream.clear();
        _lineStream.str(*it);

        _lineStream >> iatom >> fatom >> agroup;
        if (_lineStream.fail()){
        io::messages.add("bad line in GAMDATOMS block; expected format:\nGAMDATOMS\n#  INATOM   FINATOM   AGROUP\nEND\n",
                            "In_GAMD", io::message::error);
        return;
        }

        // convert to gromos
        --iatom; --fatom;

        if (iatom > topo.num_atoms() || fatom > topo.num_atoms()){
            std::ostringstream msg;
            msg << "GAMD atoms (" << iatom+1 << "-" << fatom+1 
            << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_GAMD", io::message::error);
            return;
        }

        if (agroup >= sim.param().gamd.agroups){
            std::ostringstream msg;
            msg << "GAMD atoms acceleration group out of range.";
            io::messages.add(msg.str(), "In_GAMD", io::message::error);
            return;

        }
        //load atoms
        for (unsigned int i = iatom; i <= fatom; i++){
            topo.gamd_accel_group()[i] = agroup;
        }
    } // loop over lines
  } // GAMD ATOMS
}