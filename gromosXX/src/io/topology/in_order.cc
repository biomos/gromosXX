/**
 * @file in_order.cc
 * implements methods of In_Orderparamresspec
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_order.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section orderparamresspec ORDERPARAMRESSPEC block
 * The ORDERPARAMRESSPEC block is read from the order parameter restraints specification
 * file.
 *
 * \c DISH is the carbon-hydrogen, \c DISC the carbon-carbon distance.
 * See @ref util::virtual_type for valid virtual atom types.
 * \c RN is the normalisation distance
 * \c S0 is the experimental order parameter
 * \c DS0 is a deviation which is allowed
 * \c WOPR a weight factor (multiplied by the force  constant specified in the input file
 * @verbatim
ORDERPARAMRESSPEC
# DISH  DISC
  0.1   0.153
# i  j  k  l  type    i  j  k  l  type    RN     S0    DS0   WOPR
  1  0  0  0  0       10 12 11 0   3       0.1    0.8   0.1   1.0
END
@endverbatim
 * @sa util::virtual_type util::Virtual_Atom
 */
void
io::In_Orderparamresspec::read(topology::Topology& topo,
        simulation::Simulation & sim,
        std::ostream & os) {
  std::set<std::string> block_read;
  DEBUG(7, "reading in a order parameter restraints file");

  std::vector<std::string> buffer;

  { // ORDERPARAMRESSPEC
    DEBUG(10, "ORDERPARAMRESSPEC block");
    buffer = m_block["ORDERPARAMRESSPEC"];
    block_read.insert("ORDERPARAMRESSPEC");
    if (buffer.size() <= 2) {
      io::messages.add("no or empty ORDERPARAMRESSPEC block in order param restraints file",
              "In_Orderparamresspec", io::message::error);
      return;
    } else {
      std::vector<std::string>::const_iterator it = buffer.begin() + 1,
              to = buffer.end() - 1;

      double dish = 0.0, disc = 0.0;
      bool nr_atoms = true;

      DEBUG(10, "reading in ORDERPARAMRESSPEC data");

      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> dish >> disc;

      if (_lineStream.fail()) {
        std::ostringstream msg;
        msg << "bad line in ORDERPARAMRESSPEC block: failed to read in DISH and DISC" << std::endl;
        io::messages.add(msg.str(), "In_Orderparamresspec", io::message::error);
      }

      ++it;
      for (unsigned int line_number = 2; it != to; ++line_number, ++it) {
        DEBUG(11, "\tnr " << line_number - 2);
        int type1 = 0, type2 = 0;
        std::vector<int> atom1, atom2;
        double dist_norm = 0.0, S0 = 0.0, dS0 = 0.0, w = 0.0;

        _lineStream.clear();
        _lineStream.str(*it);

        for (unsigned int i = 0; i < io::In_Orderparamresspec::MAX_ATOMS; i++) {
          int atom = 0;
          _lineStream >> atom;
          // -1 because we directly convert to array indices
          if (atom > 0) {
            atom1.push_back(atom - 1);
          } else if (atom < 0) {
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than "
                    << io::In_Orderparamresspec::MAX_ATOMS << " atoms" << std::endl;
            io::messages.add(msg.str(), "In_Orderparamresspec", io::message::error);
            nr_atoms = false;
          }
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Orderparamresspec::MAX_ATOMS; i++) {
          int atom = 0;
          _lineStream >> atom;
          if (atom > 0) {
            atom2.push_back(atom - 1);
          } else if (atom < 0) {
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than "
                    << io::In_Orderparamresspec::MAX_ATOMS << " atoms" << std::endl;
            io::messages.add(msg.str(), "In_Orderparamresspec", io::message::error);
            nr_atoms = false;
          }
        }
        _lineStream >> type2;
        _lineStream >> dist_norm >> S0 >> dS0 >> w;

        if (_lineStream.fail()) {
          std::ostringstream msg;
          msg << "bad line in ORDERPARAMRESSPEC block: " << line_number << std::endl
                  << "          " << *it;
          io::messages.add(msg.str(), "In_Orderparamresspec", io::message::error);
        }
        
        if (dist_norm <= 0.0) {
          io::messages.add("ORDERPARAMRESSPEC block: RN has to be > 0.0",
                  "In_Orderparamresspec", io::message::error);
        }
        
        if (S0 <= 0.0) {
          io::messages.add("ORDERPARAMRESSPEC block: S0 has to be > 0.0",
                  "In_Orderparamresspec", io::message::error);
        }
        
        if (dS0 <= 0.0) {
          io::messages.add("ORDERPARAMRESSPEC block: DS0 has to be > 0.0",
                  "In_Orderparamresspec", io::message::error);
        }
        
        if (w < 0.0) {
          io::messages.add("ORDERPARAMRESSPEC block: W has to be >= 0.0",
                  "In_Orderparamresspec", io::message::error);
        }

        if (nr_atoms) {
          util::virtual_type t1 = util::virtual_type(type1);
          util::virtual_type t2 = util::virtual_type(type2);

          util::Virtual_Atom v1(t1, atom1, dish, disc);
          util::Virtual_Atom v2(t2, atom2, dish, disc);


          topo.order_parameter_restraints().push_back
                  (topology::order_parameter_restraint_struct(v1, v2, dist_norm, S0, dS0, w));
        }
      } // for lines
    } // if block
  } // ORDERPARAMRESSPEC


  for (std::map<std::string, std::vector<std::string> >::const_iterator
    it = m_block.begin(),
          to = m_block.end(); it != to; ++it) {

    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add("block " + it->first + " not supported!",
              "In_Orderparamresspec", io::message::error);
    }
  }
}

