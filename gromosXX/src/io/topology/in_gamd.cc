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
 * The GAMDATOMS specifies the atoms that define each of the acceleration groups.
 * Atoms cannot be assigned to more than one group.
 * All non specified atoms will be considered to belong to group 0
 * 
 * The GAMDGROUPS defines how each of the defined groups interact with each other and with themselves.
 * The acceleration group 0 is the non accelerated one.
 * Any combination of groups that it is not specified will be considered not to be accelerated, (ACCELGROUP = 0)  
 *
 * Both blocks need to have specified how many non zero different groups are there to allocate the arrays correctly
 * The block is read from the gamd specification file
 * (\@gamd).
 *
 * @verbatim
GAMDATOMS
3
#  INATOM   FINATOM   GROUP
   1        230       0
   300      400       1
   500      800       2
GAMDGROUPS
3
# GROUP_1   GROUP_2   ACCELGROUP
  0         0         0
  1         0         1
  2         0         2
  1         1         1
  1         2         3
  2         2         2   
END
@endverbatim
 */
void
io::In_GAMD::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a GAMD definition file");

  std::vector<std::string> buffer;
  std::vector<unsigned int> atracker;

  { // GAMDATOMS

    buffer = m_block["GAMDATOMS"];
    DEBUG(10, "GAMDATOMS block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no GAMDATOMS block in GAMD definition file",
		       "In_GAMD", io::message::error);
      return;
    }

    unsigned int numgroups;
    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> numgroups;

    if (_lineStream.fail()){
        io::messages.add("bad line in GAMDATOMS block; expected atom GROUPS\n",
                            "In_GAMD", io::message::error);
        return;
    }

    if (sim.param().gamd.agroups != numgroups){
        io::messages.add("Inconsistent number of groups between imd file and gamd file",
                          "In_GAMD", io::message::error);
        return;
    }
    // check parameters to be sure numgroups have equal size in the gamd file and the imd
    ++it;

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
        DEBUG(1, "gamd groups" << numgroups);
        if (agroup > numgroups){
            std::ostringstream msg;
            msg << "GAMD atoms acceleration group out of range.";
            io::messages.add(msg.str(), "In_GAMD", io::message::error);
            return;

        }
        //load atoms
        for (unsigned int i = iatom; i <= fatom; i++){
          if(std::find(atracker.begin(), atracker.end(), i) != atracker.end()) {
            std::ostringstream msg;
            msg << "Duplicated atoms in the GAMDATOMS.";
            io::messages.add(msg.str(), "In_GAMD", io::message::error);
            return;
            } else {
              topo.gamd_accel_group()[i] = agroup;
              atracker.push_back(i);
            }
        }
    } // loop over lines
  } // GAMD ATOMS

  { // GAMDGROUPS
    buffer = m_block["GAMDGROUPS"];
    DEBUG(10, "GAMDATOMS block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no GAMDGROUPS block in GAMD definition file",
		       "In_GAMD", io::message::error);
      return;
    }
    unsigned int anumgroups;
    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> anumgroups;

    if (_lineStream.fail()){
        io::messages.add("bad line in GAMDATOMS block; expected acceleration GROUPS\n",
                            "In_GAMD", io::message::error);
        return;
    }

    if (sim.param().gamd.igroups != (anumgroups + 1)){
        io::messages.add("Inconsistent number of groups between imd file and gamd file",
                          "In_GAMD", io::message::error);
        return;
    }
    // check parameters to be sure numgroups have equal size in the gamd file and the imd
    ++it;

    DEBUG(10, "reading in GAMDGROUPS data");
    unsigned int n, group1, group2, intgroup;
    for(n=0; it != to; ++n, ++it){
      DEBUG(11, "\tnr " << n);

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> group1 >> group2 >> intgroup;
      if (_lineStream.fail()){
      io::messages.add("bad line in GAMDGROUPS block; expected format:\nGAMDGROUPS\n#  GROUP_1   GROUP_2   IGROUP\nEND\n",
                          "In_GAMD", io::message::error);
      return;
      }
      if(intgroup > anumgroups){
      io::messages.add("Acceleration groups in GAMDGROUPS are out of range",
                          "In_GAMD", io::message::error);
      return;        
      }
      // Adds both possible combinations
      std::vector<unsigned int> newkey1 = {group1, group2};
      std::vector<unsigned int> newkey2 = {group2, group1};
      topo.get_gamd_interaction_pairs().insert({newkey1, intgroup});
      topo.get_gamd_interaction_pairs().insert({newkey2, intgroup});
    } // loop over lines  
  } // GAMDGROUPS
}
