/**
 * @file in_zaxisoribias.cc
 * implements methods of In_Topology.
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_zaxisoribias.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;


void
io::In_Zaxisoribias::read(topology::Topology& topo,
              simulation::Simulation & sim,
              std::ostream & os){

  DEBUG(7, "reading in a z-axis orientation bias file");

  if (!quiet)
    os << "Z-AXIS ORIENTATION BIAS\n";
  read_ZAXISORIBIASSPEC(topo, sim, os);

  if (!quiet) os << "END\n";

  if (sim.param().zaxisoribias.zaxisoribias) {
    if (block_read.count("ZAXISORIBIASSPEC") == 0)
      io::messages.add("Z-axis orientation bias on but no ZAXISORIBIASSPEC found",
               "In_Zaxisoribias",
               io::message::error);
  }

  for(std::map<std::string, std::vector<std::string> >::const_iterator
    it = m_block.begin(),
    to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add(" block " + it->first + " not supported!",
               "In_Zaxisoribias",
               io::message::warning);
    }
  }

}

/**
  * @section zaxisoribiasspec ZAXISORIBIASSPEC block
  * This block is for biasing the angle of a vector with the z-axis. It is read
  * from the z-axis bias specifcation file.
 * @snippet snippets/snippets.cc ZAXISORIBIASSPEC
 */

void io::In_Zaxisoribias::read_ZAXISORIBIASSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os)  { // DISTANCERES
  DEBUG(10, "ZAXISORIBIASSPEC block");

  std::stringstream exampleblock;

  //###BARTOSZ### What has to be changed here?

  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "ZAXISORIBIASSPEC\n";
  exampleblock << "# DISH, DISC carbon-hydrogen/carbon-carbon distance\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use)\n";
  exampleblock << "# type   virtual atom type\n";
  exampleblock << "# a0, w0  target z-axis angle and force constant weighting factor\n";
  exampleblock << "# rah    form and dimension of the potential\n";
  exampleblock << "# full harmonic:\n";
  exampleblock << "#     0: align along z\n";
  exampleblock << "#        (no other values implemented at the moment)\n";
  exampleblock << "#  repulsive or attractive potential\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# i  j  k  l  type    i  j  k  l  type    a0    w0    rah\n";
  exampleblock << "  1  0  0  0  0       10 12 11 13 3       45   1.0   0\n";
  exampleblock << "END\n";

  //###BARTOSZ### What has to be changed here? ^

  std::string blockname = "ZAXISORIBIASSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    double dish,disc;

    DEBUG(10, "reading in ZAXISORIBIAS data");

    if (!quiet) {
        switch (sim.param().zaxisoribias.zaxisoribias) {
          case 0:
            os << "\tZ-axis orientation bias OFF\n";
            // how did you get here?
            break;
          case 1:
            os << "\tZ-axis orientation bias ON\n";
            break;
          case -1:
            os << "\tZ-axis orientation bias ON\n"
                    << "\tno arccos function involved\n";
            break;
          case 2:
            os << "\tZ-axis orientation bias ON\n"
                    << "\t\t(using force constant K*w0)\n";
            break;
          case -2:
            os << "\tZ-axis orientation bias ON\n"
                    << "\tno trigonometric function involved\n";
            break;
          default:
            os << "\tZ-axis orientation bias ERROR\n";
        }
    }

    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

    if (!quiet){
    os << std::setw(10) << "DISH"
       << std::setw(10) << "DISC"
       << "\n"
       <<  std::setw(10)<< dish
       <<  std::setw(10)<< disc
       << "\n";

    os << std::setw(10) << "i"
       << std::setw(8) << "j"
       << std::setw(8) << "k"
       << std::setw(8) << "l"
       << std::setw(5) << "type"
       << std::setw(10) << "i"
       << std::setw(8) << "j"
       << std::setw(8) << "k"
       << std::setw(8) << "l"
       << std::setw(5) << "type"
       << std::setw(8) << "a0"
       << std::setw(8) << "w0"
       << std::setw(4) << "rah"
       << "\n";
      }

    unsigned int num = block.numlines()-3;
    for(unsigned int line_number=0; line_number < num; line_number++){

      DEBUG(11, "\tnr " << line_number);

      int type1, type2;
      std::vector<int> atom1, atom2;
      double a0,w0;
      int rah;

      for(unsigned int i = 0; i < io::In_Zaxisoribias::MAX_ATOMS; i++) {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM["+str_i+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Zaxisoribias", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for(unsigned int i = 0; i < io::In_Zaxisoribias::MAX_ATOMS; i++) {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM["+str_i+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Zaxisoribias", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      block.get_next_parameter("A0", a0, ">=0", "");
      block.get_next_parameter("W0", w0, ">=0", "");
      block.get_next_parameter("RAH", rah, "", "");

      // move into range 0 to 180 degrees
        while (a0 > 360) a0 -= 360;
        if (a0 < 0) a0 *= -1;

    // g++ 3.2 fix
      if(!block.error()){
        util::virtual_type t1 = util::virtual_type(type1);
        util::virtual_type t2 = util::virtual_type(type2);

        util::Virtual_Atom v1(t1, atom1, dish, disc);
        util::Virtual_Atom v2(t2, atom2, dish, disc);

        topo.zaxisori_restraints().push_back
          (topology::zaxisori_restraint_struct(v1,v2,a0 * math::Pi / 180,w0,rah));

        if (!quiet){
          for(unsigned int i = 0; i < io::In_Zaxisoribias::MAX_ATOMS; i++) {
          // the first element has the width 10, if i is bigger then the number of atoms
          // specified, just print 0.
            os << std::setw(i == 0 ? 10 : 8) << (i < atom1.size() ? atom1[i]+1 : 0);
          }
          os << std::setw(5) << type1;

          for(unsigned int i = 0; i < io::In_Zaxisoribias::MAX_ATOMS; i++) {
            os << std::setw(i == 0 ? 10 : 8) << (i < atom2.size() ? atom2[i]+1 : 0);
          }
          os << std::setw(5) << type2
             << std::setw(8) << a0
             << std::setw(8) << w0
             << std::setw(4) << rah
             << "\n";
        }
      }
    } // for restraints
    block.get_final_messages();
  } // if block content
} // ZAXISORIBIAS