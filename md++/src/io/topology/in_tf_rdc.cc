/**
 * @file in_tf_rdc.cc
 * implements methods of In_Tfrdcresspec
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_tf_rdc.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section tfrdcresspec TFRDCRESSPEC block
 * The TFRDCRESSPEC block is read from the tensor-free RDC restraints
 *  specification file .tfr.
 *
 * \c DISH is the carbon-hydrogen, \c DISC the carbon-carbon distance.
 * See @ref util::virtual_type for valid virtual atom types.
 * \c R0 is the normalisation distance
 * \c G1/G2 is the gyromagnetic ratio for nucleus 1/nucleus 2 in
 * 10^6 rad /(T s)
 * \c D0 is the experimental RDC in Hz
 * \c DD0 is a deviation which is allowed (in Hz)
 * \c WTFRDC a weight factor (multiplied by the force constant specified in
 * the input file)
 * @verbatim
TFRDCRESSPEC
#DISH  DISC
   0.1   0.13
# i   j   k   l   type  i    j    k   l  type  R0      G1      G2      D0      DD0   WTFRDC
  2   0   0   0   0     4    0    0   0    0   0.1  267.513  19.331   0.136    0.05     1.0
  3   0   0   0   0     5    0    0   0    0   0.1  267.513  67.262   0.136    0.05     1.0
END
@endverbatim
 * @sa util::virtual_type util::Virtual_Atom
 */
void
io::In_Tfrdcresspec::read(topology::Topology& topo,
        simulation::Simulation & sim,
        std::ostream & os) {
  std::set<std::string> block_read;
  DEBUG(7, "reading in a tensor-free RDC restraints file");

  std::vector<std::string> buffer;

  { // TFRDCRESSPEC
    DEBUG(10, "TFRDCRESSPEC block");
    buffer = m_block["TFRDCRESSPEC"];
    block_read.insert("TFRDCRESSPEC");
    if (buffer.size() <= 2) {
      io::messages.add("no or empty TFRDCRESSPEC block in tensor-free "
         "RDC restraints file", "In_Tfrdcresspec", io::message::error);
      return;
    } else {
      std::vector<std::string>::const_iterator it = buffer.begin() + 1,
              to = buffer.end() - 1;

      double dish, disc;
      bool nr_atoms = true;

      DEBUG(10, "reading in TFRDCRESSPEC data");

      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> dish >> disc;

      if (_lineStream.fail()) {
        std::ostringstream msg;
        msg << "bad line in TFRDCRESSPEC block: failed to read in DISH and DISC" << std::endl;
        io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
      }

      ++it;
      for (unsigned int line_number = 2; it != to; ++line_number, ++it) {
        DEBUG(11, "\tnr " << line_number - 2);
        int type1, type2;
        std::vector<int> atom1, atom2;
        double dist_norm, gyri, gyrj, D0, dD0, w;

        _lineStream.clear();
        _lineStream.str(*it);

        for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; i++) {
          int atom;
          _lineStream >> atom;
          // -1 because we directly convert to array indices
          if (atom > 0) {
            atom1.push_back(atom - 1);
          } else if (atom < 0) {
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than "
                    << io::In_Tfrdcresspec::MAX_ATOMS << " atoms" << std::endl;
            io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
            nr_atoms = false;
          }
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; i++) {
          int atom;
          _lineStream >> atom;
          if (atom > 0) {
            atom2.push_back(atom - 1);
          } else if (atom < 0) {
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than "
                    << io::In_Tfrdcresspec::MAX_ATOMS << " atoms" << std::endl;
            io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
            nr_atoms = false;
          }
        }
        _lineStream >> type2;
        _lineStream >> dist_norm >> gyri >> gyrj >> D0 >> dD0 >> w;

        // Scaling gyromagnetic ratios from 10^6 rad/(T*s) to rad /(e*u)
        gyri *= 0.010364272;
        gyrj *= 0.010364272;
        // Scaling RDC and DeltaRDC from Hz to 1/ps
        D0 *= 0.000000000001;
        dD0 *= 0.000000000001;
        if (_lineStream.fail()) {
          std::ostringstream msg;
          msg << "bad line in TFRDCRESSPEC block: " << line_number << std::endl
                  << "          " << *it;
          io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
        }

        if (dist_norm <= 0.0) {
          io::messages.add("TFRDCRESSPEC block: R0 has to be > 0.0",
                  "In_Tfrdcresspec", io::message::error);
        }

        if(fabs(-( math::eps0_i * math::h_bar * gyri * gyrj )/( pow(math::spd_l,2) * 4.0 * pow(math::Pi,2) * pow(dist_norm,3.0))) < fabs(D0)){
          io::messages.add("The chosen RDC is larger in magnitude than RDC_max. This is probably a mistake and may result in strange behaviour.",
           "In_Tfrdcresspec", io::message::warning);
        }

        if (dD0 < 0.0) {
          io::messages.add("TFRDCRESSPEC block: DD0 has to be >= 0.0",
                  "In_Tfrdcresspec", io::message::error);
        }

        if (w < 0.0) {
          io::messages.add("TFRDCRESSPEC block: WTFRDC has to be >= 0.0",
                  "In_Tfrdcresspec", io::message::error);
        }

        if (nr_atoms) {
          util::virtual_type t1 = util::virtual_type(type1);
          util::virtual_type t2 = util::virtual_type(type2);

          util::Virtual_Atom v1(t1, atom1, dish, disc);
          util::Virtual_Atom v2(t2, atom2, dish, disc);


          topo.tf_rdc_restraints().push_back
                  (topology::tf_rdc_restraint_struct(v1, v2, dist_norm, gyri, gyrj, D0, dD0, w));
        }
      } // for lines
    } // if block
} // TFRDCRESSPEC


  for (std::map<std::string, std::vector<std::string> >::const_iterator
    it = m_block.begin(),
          to = m_block.end(); it != to; ++it) {

    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add("block " + it->first + " not supported!",
              "In_Tfrdcresspec", io::message::error);
    }
  }
}
