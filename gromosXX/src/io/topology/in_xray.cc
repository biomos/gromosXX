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
@endverbatim
 *
 * @section xrayrespara XRAYRESPARA block
 * The XRAYRESPARA specifies the parameters needed to calculate the structure
 * factors with clipper-library:
 *  - Spacegroup SPGR in Hermann-Mauguin-Form
 *  - Unitcell CELL in form: a b c alpha beta gamma
 *  - Scattering-Resolution RESO in nm
 *  - Standard-Bfactor STDBF used for all atoms
 *
 * The block is read from the xray restraints specification file
 * (\@xray).
 *
 * @verbatim
XRAYRESPARA
#  SPGR
   P 21 21 21
#  CELL                                            RESO    STDBF
   5.9062  6.8451  3.0517  90.00  90.00  90.00     0.15      1.0
END
@endverbatim
 */
void
io::In_Xrayresspec::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a xray restraints file");

  if (!quiet)

    os << "XRAY RESTRAINTS\n";

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

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

    unsigned int i, h, k, l;
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
    _lineStream >> sim.param().xrayrest.cell_a >> sim.param().xrayrest.cell_b
            >> sim.param().xrayrest.cell_c >> sim.param().xrayrest.cell_alpha
            >> sim.param().xrayrest.cell_beta
            >> sim.param().xrayrest.cell_gamma
            >> sim.param().xrayrest.resolution
            >> sim.param().xrayrest.bfactor ;

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

  /* OBSOLETE
  { // XRAYFITATOMS

    buffer = m_block["XRAYFITATOMS"];
    DEBUG(10, "XRAYFITATOMS block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no XRAYFITATOMS block in position restraints file",
		       "In_Xrayresspec", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in XRAYFITATOMS data");

    unsigned int i, nr;
    math::Vec pos;
    for (i = 0; it != to; ++i, ++it) {

      DEBUG(11, "\tnr " << i);

      std::string line(*it);
      if (line.length() < 17) {
        io::messages.add("line too short in XRAYFITATOMS block", "In_Xrayresspec",
                io::message::error);
      }

      // the first 17 chars are ignored
      line.erase(line.begin(), line.begin() + 17);

      _lineStream.clear();
      _lineStream.str(line);

      _lineStream >> nr >> pos(0) >> pos(1) >> pos(2);

      if (_lineStream.fail()) {
        io::messages.add("bad line in XRAYFITATOMS block",
                "In_Xrayresspec", io::message::error);
        return;
      }

      if (nr>topo.num_solute_atoms()) {
        io::messages.add("bad line in XRAYFITATOMS block: atom number out of range",
                "In_Xrayresspec", io::message::error);
        return;
      }

      topo.xray_fitatoms().push_back(topology::xray_fitatom_struct(nr-1, pos));

    }
  } // XRAYFITATOMS
  */


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

      case simulation::xrayrest_loel :
                os << "\tXray local-elevation restraints ON\n";
        break;
    }

    os << "END\n";
  }
}


