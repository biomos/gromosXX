/**
 * @file in_xray.cc
 * implements methods of In_xrayresspec and In_xrayres.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>
#include <io/configuration/in_configuration.h>
#include <io/configuration/out_configuration.h>

#include "in_xray.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section xrayresspec XRAYRESSPEC block
 * The XRAYRESSPEC specifies the experimental structure factors on which we
 * restrain.
 *
 * The block is read from the xray restraints specification file
 * (\@xray).
 *
 * @verbatim
XRAYRESSPEC
#  H   K   L          SF   STDDEV_SF
   0   0   2       13.46        0.22
   0   0   4      184.62        1.34
   0   0   6       20.99        0.33
   0   0   8      158.82        1.16
   0   0  10       80.81        0.64
END
@endverbatim
 *
 * @section elementmap XRAYELEMENTSPEC block
 * The XRAYELEMENTSPEC block contains the elementname of all atoms.
 *
 * The block is read from the xray restraints specification file
 * (\@xray).
 *
 * @verbatim
XRAYELEMENTSPEC
#  ELEMENT[1..N]
N N N H H O N C Na Ca
H O N C Na Ca H O N C
END
XRAYSOLVELEMENTSPEC
#  ELEMENT[1..N]
O H H
END
@endverbatim
 *
 * @section xraybfoccspec XRAYBFOCCSPEC block
 * The XRAYBFOCCSPEC block contains the B-factors and occupancies of all atoms.
 *
 * The block is read from the xray restraints specification file (\@xray).
 *
 * @verbatim
XRAYBFOCCSPEC
#  BF[1..N]  OCC[1..N]
   0.01      1.0
   0.02      0.9
END
XRAYSOLVBFOCCSPEC
#  BF[1..N]  OCC[1..N]
   0.01      1.0
   0.01      0.0
   0.01      0.0
END
@endverbatim
 *
 * @section xrayrespara XRAYRESPARA block
 * The XRAYRESPARA specifies the parameters needed to calculate the structure
 * factors with clipper-library:
 *  - Spacegroup SPGR in Hermann-Mauguin-Form
 *  - Scattering-Resolution RESO in nm
 *
 * The block is read from the xray restraints specification file
 * (\@xray).
 *
 * @verbatim
XRAYRESPARA
#  SPGR
   P 21 21 21
#  RESO
   0.15
END
@endverbatim
 *
 * @section xrayumbrellaweight XRAYUMBRELLAWEIGHT block
 * The XRAYUMBRELLAWEIGHT specified the atoms that are attached to a particular
 * umbrella. The umbrella is than weighted by the electron density deviation instead
 * of the number of visits.
 *
 * The block is read from the xray restraints specification file (\@xray)
 * @verbatim
XRAYUMBRELLAWEIGHT
# UMBID      : ID of the umbrella
# THRES      : Threshold for the flat bottom potatial
# DENSCUT    : Cutoff for integration over atom volume
# ATOMS[1..N]: The atoms used for weighting the umbrella
1 0.23 0.5 23 24 24 26
2 0.33 0.5 56 57 58 59 60
END
@endverbatim
 *
 */
void
io::In_Xrayresspec::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a xray restraints file");

  if (!quiet)

    os << "XRAY RESTRAINTS\n";

  std::vector<std::string> buffer;

  { // XRAYRESSPEC

    buffer = m_block["XRAYRESSPEC"];
    DEBUG(10, "XRAYRESSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYRESSPEC block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in XRAYRESSPEC data");

    int i, h, k, l;
    double sf, stddev_sf;
    for(i=0; it != to; ++i, ++it){

      DEBUG(11, "\tnr " << i);

      std::string line(*it);

      _lineStream.clear();
      _lineStream.str(line);

      _lineStream >> h >> k >> l >> sf >> stddev_sf;

      DEBUG(11, "grid indices " << h << " " << k << " " << l << " St_Fac's "
              << sf << " " << stddev_sf);

      if (_lineStream.fail()){
        io::messages.add("bad line in XRAYRESSPEC block",
                "In_Xrayresspec", io::message::error);
        return;
      }

      topo.xray_restraints().push_back(topology::xray_restraint_struct(h, k, l, sf, stddev_sf));

    }
  } // XRAYRESSPEC

  { // XRAYRESPARA

    buffer = m_block["XRAYRESPARA"];
    DEBUG(10, "XRAYRESPARA block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYRESPARA block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }

    if (buffer.size()!=4){
      io::messages.add("XRAYRESPARA block must contain 2 lines",
		       "in_Xrayresspec", io::message::error);
      return;
    }

    DEBUG(10, "reading in XRAYRESPARA data");

    std::string spacegroup = buffer[1];
    const size_t startpos = spacegroup.find_first_not_of(" \t"); // start trim
    const size_t endpos = spacegroup.find_last_not_of(" \t"); // end trim
    spacegroup = spacegroup.substr( startpos, endpos-startpos+1 );
    sim.param().xrayrest.spacegroup=spacegroup;
    // TODO: Check wether valid or not (-> clipper_message.h)

    _lineStream.clear();
    _lineStream.str(buffer[2]);
    _lineStream >> sim.param().xrayrest.resolution;

    if (_lineStream.fail()) {
      io::messages.add("bad second line in XRAYRESPARA block",
              "In_Xrayresspec", io::message::error);
      return;
    }
  } // XRAYRESPARA


  { // XRAYELEMENTSPEC

    buffer = m_block["XRAYELEMENTSPEC"];
    DEBUG(10, "XRAYELEMENTSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYELEMENTSPEC block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    topo.xray_elements().resize(topo.num_solute_atoms(), "H");
    for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {;
      _lineStream >> topo.xray_elements()[i];
      
      if (_lineStream.fail()) {
        topo.xray_elements()[i] = "H";
        io::messages.add("bad line in XRAYELEMENTSPEC block",
                "In_Xrayresspec", io::message::error);
        return;
      }
    }

  } // XRAYELEMENTSPEC

  { // XRAYSOLVELEMENTSPEC

    buffer = m_block["XRAYSOLVELEMENTSPEC"];
    DEBUG(10, "XRAYSOLVELEMENTSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYSOLVELEMENTSPEC block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    topo.xray_solvelements().resize(topo.solvent(0).num_atoms(), "H");
    for (unsigned int i = 0; i < topo.solvent(0).num_atoms(); ++i) {;
      _lineStream >> topo.xray_solvelements()[i];

      if (_lineStream.fail()) {
        topo.xray_solvelements()[i] = "H";
        io::messages.add("bad line in XRAYELEMENTSPEC block",
                "In_Xrayresspec", io::message::error);
        return;
      }
    }
  } // XRAYSOLVELEMENTSPEC

  { // XRAYBFOCCSPEC
    buffer = m_block["XRAYBFOCCSPEC"];
    DEBUG(10, "XRAYBFOCCSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYBFOCCSPEC block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    topo.xray_b_factors().resize(topo.num_solute_atoms(), 0.01);
    topo.xray_occupancies().resize(topo.num_solute_atoms(), 0.0);
    for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {;
      _lineStream >> topo.xray_b_factors()[i] >> topo.xray_occupancies()[i];

      if (_lineStream.fail()) {
        topo.xray_b_factors()[i] = 0.01;
        topo.xray_occupancies()[i] = 0.0;
        io::messages.add("bad line in XRAYBFOCCSPEC block",
                "In_Xrayresspec", io::message::error);
        return;
      }
    }
  } // XRAYBFOCCSPEC

  { // XRAYSOLVBFOCCSPEC
    buffer = m_block["XRAYSOLVBFOCCSPEC"];
    DEBUG(10, "XRAYSOLVBFOCCSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYSOLVBFOCCSPEC block in xray restraints file",
		       "in_Xrayresspec", io::message::error);
      return;
    }
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    topo.xray_solv_b_factors().resize(topo.solvent(0).num_atoms(), 0.01);
    topo.xray_solv_occupancies().resize(topo.solvent(0).num_atoms(), 0.0);
    for (unsigned int i = 0; i < topo.solvent(0).num_atoms(); ++i) {;
      _lineStream >> topo.xray_solv_b_factors()[i] >> topo.xray_solv_occupancies()[i];

      if (_lineStream.fail()) {
        topo.xray_solv_b_factors()[i] = 0.01;
        topo.xray_solv_occupancies()[i] = 0.0;
        io::messages.add("bad line in XRAYSOLVBFOCCSPEC block",
                "In_Xrayresspec", io::message::error);
        return;
      }
    }
  } // XRAYSOLVBFOCCSPEC

  { // XRAYUMBRELLAWEIGHT
    buffer = m_block["XRAYUMBRELLAWEIGHT"];
    DEBUG(10, "XRAYUMBRELLAWEIGHT block : " << buffer.size());

    if (sim.param().xrayrest.local_elevation) {
      if (!buffer.size()) {
        io::messages.add("no XRAYUMBRELLAWEIGHT block in xray restraints file",
                "in_Xrayresspec", io::message::error);
        return;
      }
      std::vector<std::string>::const_iterator it = buffer.begin() + 1,
              to = buffer.end() - 1;

      for(; it != to; ++it) {
        _lineStream.clear();
        // trim whitespace from right end
        std::string line(*it);
        std::string::size_type right = line.find_last_not_of(" \n\r\t");
        if (right != std::string::npos)
          line.erase(right+1);
        _lineStream.str(line);
        int umbrella;
        double thres, cut;
        _lineStream >> umbrella >> thres >> cut;
        if (_lineStream.fail()) {
          io::messages.add("XRAYUMBRELLAWEIGHT block: Cannot read umbrella id and threshold",
                  "in_Xrayresspec", io::message::error);
          return;
        }
        if (thres < 0.0) {
          io::messages.add("XRAYUMBRELLAWEIGHT block: Bad threshold (<0.0)",
                  "in_Xrayresspec", io::message::error);
          return;
        }
        std::vector<unsigned int> atoms;
        int atom;
        while(!_lineStream.eof()) {
          _lineStream >> atom;
          --atom;
          if (_lineStream.fail() || atom < 0 || unsigned(atom) > topo.num_atoms()) {
            io::messages.add("XRAYUMBRELLAWEIGHT block: Bad atom number",
                  "in_Xrayresspec", io::message::error);
            return;
          }
          atoms.push_back(atom);
        } // while atoms
        if (atoms.empty()) {
          io::messages.add("XRAYUMBRELLAWEIGHT block: no atoms given",
                  "in_Xrayresspec", io::message::error);
          return;
        }

        topo.xray_umbrella_weights().push_back(topology::xray_umbrella_weight_struct(umbrella, thres, cut, atoms));

      } //  for lines
    } else {
      if (buffer.size()) {
        io::messages.add("XRAYUMBRELLAWEIGHT block in xray restraints file not read in",
                "in_Xrayresspec", io::message::warning);
      }
    }
  } // XRAYSOLVBFOCCSPEC

  if (!quiet) {
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
                os << "\tXray restraints OFF\n";
        break;
      case simulation::xrayrest_inst :
                os << "\tXray instantaneous restraints ON\n";
        break;
      case simulation::xrayrest_avg :
                os << "\tXray time-averaged restraints ON\n";
        break;
      case simulation::xrayrest_biq :
                os << "\tXray biquadratic instantaneous/time-averaged restraints ON\n";
        break;
    }
    os << "END\n";
  }
}


