/**
 * @file in_xray.cc
 * implements methods of In_xrayresspec and In_xrayres.
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

#include "in_xray.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section xrayresspec XRAYRESSPEC block
 * The XRAYRESSPEC and XRAYRFREESPEC block specify the experimental structure factors on which we
 * restrain or use the calculate R-free
 *
 * The blocks are read from the xray restraints specification file
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
XRAYRFREESPEC
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
 *  - Conversion factor for the length unit to Angstrom (TOANG)
 *
 * The block is read from the xray restraints specification file
 * (\@xray).
 *
 * @verbatim
XRAYRESPARA
#  SPGR
   P 21 21 21
#  RESO  TOANG
   0.15   10.0
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
 * @section xraysymresspec XRAYSYMRESSPEC block
 * The XRAYSYMRESSPEC block specifies parameters of the method. On two additional
 * lines data about the spacegroup and asymmetric units have to be given. The rest is
 * just specification of the atoms.
 *
 * @verbatim
XRAYSYMRESSPEC
# NTSYM      : use symmetry restraints? (default 0)
#              - 0: do not use symmetry restraints
#              - 1: restrain ASUs in a pairwise way
#              - 2: constrain atoms of additional ASUs to image of first ASU
# CSYM       : force constant for symmetry restraints (default 0)
# SYMSPGR    : spacegroup for symmetry restraints (default 0)
# ASUDEF     : pointer to first atom in every ASU
# SYMATOMS   : the atoms to restrain
#
# NTSYM     CSYM
      1     25.0
# SYMSPGR
  P 21 21 2
# ASUDEF[1..SYMNUMSYM]
  1  101  201  301
# SYMATOMS
    2 HEXA  CH23       9
    2 HEXA  CH24      10
    2 HEXA  CH25      11
    2 HEXA  CH36      12
END
@endverbatim
 *
 * @section xraybfopt XRAYBFACTOROPTIMISATION block
 * The XRAYBFACTOROPTIMISATION block is used to specify the settings for
 * B factor optimisation. The optimisation is carried out by conjugate gradient
 * minimisation of the SF restraint residual.
 *
 * Groups can be defined for which the B factor is identical. If no groups
 * are given every atom is but into an individual group.
 *
 * @verbatim
XRAYBFACTOROPTIMISATION
# BFOPTS     : Optimise B-factors every BFOPTSth step
# BFOPTTI    : Terminate after BFOPTTI iterations
# BFOPTTG    : Terminate if gradient is smaller than BFOPTTG
# BFOPTMN    : Minimum B-factor
# BFOPTMX    : Maximum B-factor
# BFOPTNG    : The number of B factor groups.
# BFOPTGS    : Size of a group
# BFOPTGM[]  : Members of the group
#
# BFOPTS  BFOPTTI BFOPTTG BFOPTMN BFOPTMX
     100      100    0.01   0.001     1.0
# BFOPTNG
        2
# BFOPTGS  BFOPTGM[1] ...
        4    1    11    21    31
        2    2     3
END
@endverbatim
 *
 * @section xraysfcalc XRAYSFCALC block
 * The XRAYSFCALC block is used to tune the structure factor computation.
 * 
 * @verbatim
XRAYSFCALC
# SFCTOL: recalculate structure factors if an atoms has moved by SFCTOL
# SFCST: recalculate structure factors every SFCSTth step
#
# SFCTOL   SFCST
     0.1       5
END
@endverbatim
 *
 * @section xrayreplicaexchange XRAYREPLICAEXCHANGE block
 * Make resolution or force constant lambda dependent
 * @verbatim
XRAYREPLICAEXCHANGE
# NTXRRE: use X-ray replica exchange 0..2
#                   0 don't use replica exchange
#                   1 replica exchange on force constant
#                   2 replica exchange on resolution
# CXRREMN: minimal force constant / resolution
# CXRREMX: maximal force constant / resolution
#
# NTXRRE  CXREEMN CXREEMX
     1        0.0  1000.0
END
@endverbatim
 * 
 * @section xrayoverallbfactor XRAYOVERALLBFACTOR block
 * Controls the usage of an overall B-factor
 * @verbatim
XRAYOVERALLBFACTOR
# XROB >= 0.0 the overall B factor used.
# XROBF 0,1 fit the overall B factor using least-squares
#
# XROB  XROBF
          0.5            1
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

    int i = 0, h = 0, k = 0, l = 0;
    double sf = 0.0, stddev_sf = 0.0;
    topo.xray_restraints().clear();
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

      if (sf < 0.0 || stddev_sf < 0.0) {
        io::messages.add("bad line in XRAYRESSPEC block: negative SF/STDDEV",
                "In_Xrayresspec", io::message::error);
        return;
      }

      topo.xray_restraints().push_back(topology::xray_restraint_struct(h, k, l, sf, stddev_sf));

    }
  } // XRAYRESSPEC
  { // XRAYRFREESPEC
    buffer = m_block["XRAYRFREESPEC"];
    DEBUG(10, "XRAYRFREESPEC block : " << buffer.size());

    if (buffer.size()) {
      std::vector<std::string>::const_iterator it = buffer.begin() + 1,
              to = buffer.end() - 1;

      DEBUG(10, "reading in XRAYRFREESPEC data");

      int i = 0, h = 0, k = 0, l = 0;
      double sf = 0.0, stddev_sf = 0.0;
      topo.xray_rfree().clear();
      for (i = 0; it != to; ++i, ++it) {

        DEBUG(11, "\tnr " << i);

        std::string line(*it);

        _lineStream.clear();
        _lineStream.str(line);

        _lineStream >> h >> k >> l >> sf >> stddev_sf;

        DEBUG(11, "grid indices " << h << " " << k << " " << l << " St_Fac's "
                << sf << " " << stddev_sf);

        if (_lineStream.fail()) {
          io::messages.add("bad line in XRAYRFREESPEC block",
                  "In_Xrayresspec", io::message::error);
          return;
        }

        if (sf < 0.0 || stddev_sf < 0.0) {
          io::messages.add("bad line in XRAYRESSPEC block: negative SF/STDDEV",
                  "In_Xrayresspec", io::message::error);
          return;
        }

        topo.xray_rfree().push_back(topology::xray_restraint_struct(h, k, l, sf, stddev_sf));
      }
    }
  } // XRAYRFREESPEC
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

    double to_ang = 0.0;
    _lineStream >> to_ang;
    if (_lineStream.fail()) {
      io::messages.add("XRAYRESPARA block: No conversion factor for the length unit to Angstrom given. Asuming 10.0",
              "In_Xrayresspec", io::message::warning);
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
      //std::cerr << topo.xray_occupancies().size() << " xray size "<< std::endl;

   //   std::cout << " i input " << i << std::endl;
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

      topo.xray_umbrella_weights().clear();
      for(; it != to; ++it) {
        _lineStream.clear();
        // trim whitespace from right end
        std::string line(*it);
        std::string::size_type right = line.find_last_not_of(" \n\r\t");
        if (right != std::string::npos)
          line.erase(right+1);
        _lineStream.str(line);
        int umbrella = 0;
        double thres = 0.0, thres_growth = 0.0, thres_overshoot = 0.0, cut = 0.0;
        bool thres_freeze = 0;
        _lineStream >> umbrella 
                    >> thres 
                    >> thres_growth 
                    >> thres_overshoot
                    >> thres_freeze 
                    >> cut;
        if (_lineStream.fail()) {
          io::messages.add("XRAYUMBRELLAWEIGHT block: Cannot read umbrella id and threshold data",
                  "in_Xrayresspec", io::message::error);
          return;
        }
        if (thres < 0.0) {
          io::messages.add("XRAYUMBRELLAWEIGHT block: Bad threshold (<0.0)",
                  "in_Xrayresspec", io::message::error);
          return;
        }
        std::vector<unsigned int> atoms;
        int atom = 0;
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

        topo.xray_umbrella_weights().push_back(topology::xray_umbrella_weight_struct(umbrella, thres, thres_growth, thres_overshoot, thres_freeze, cut, atoms));

      } //  for lines
    } else {
      if (buffer.size()) {
        io::messages.add("XRAYUMBRELLAWEIGHT block in xray restraints file not read in",
                "in_Xrayresspec", io::message::warning);
      }
    }
  } // XRAYSOLVBFOCCSPEC
  { // XRAYSYMRESSPEC
    buffer = m_block["XRAYSYMRESSPEC"];
    DEBUG(10, "XRAYSYMRESSPEC block : " << buffer.size());

    if (!buffer.empty()) {
      if (buffer.size() < 6) {
        io::messages.add("XRAYSYMRESSPEC block: not enough lines given",
                "in_Xrayresspec", io::message::error);
        return;
      }

      _lineStream.clear();
      _lineStream.str(buffer[1]);
      unsigned int ntsym = 0;
      _lineStream >> ntsym >> sim.param().xrayrest.sym_force_constant;
      if (_lineStream.fail()) {
        io::messages.add("XRAYSYMRESSPEC block: Cannot read method and force constant from first line",
                "in_Xrayresspec", io::message::error);
        sim.param().xrayrest.sym_force_constant = 0.0;
        return;
      }
      switch (ntsym) {
        case 0:
          sim.param().xrayrest.symrest = simulation::xray_symrest_off;
          break;
        case 1:
          sim.param().xrayrest.symrest = simulation::xray_symrest_ind;
          break;
        case 2:
          sim.param().xrayrest.symrest = simulation::xray_symrest_constr;
          break;
        default:
          sim.param().xrayrest.symrest = simulation::xray_symrest_off;
          io::messages.add("XRAYSYMRESSPEC block: Invalid NTSYM",
                  "in_Xrayresspec", io::message::error);
          sim.param().xrayrest.sym_force_constant = 0.0;
          return;
      }

      std::string spacegroup = buffer[2];
      const size_t startpos = spacegroup.find_first_not_of(" \t"); // start trim
      const size_t endpos = spacegroup.find_last_not_of(" \t"); // end trim
      spacegroup = spacegroup.substr(startpos, endpos - startpos + 1);
      sim.param().xrayrest.sym_spacegroup = spacegroup;
#ifdef HAVE_CLIPPER
      clipper::Spacegroup spacegr;
      try {
        clipper::Spgr_descr spgrinit(clipper::String(sim.param().xrayrest.sym_spacegroup), clipper::Spgr_descr::HM);
        spacegr.init(spgrinit);
      } catch (const clipper::Message_fatal & msg) {
        io::messages.add("In_Xrayresspec", "XRAYSYMRESSPEC block: "+msg.text(), io::message::error);
        return;
      }
      _lineStream.clear();
      _lineStream.str(buffer[3]);
      DEBUG(8, "\tSpacegroup has " << spacegr.num_symops() << " symmetry operations.");
      if (spacegr.num_symops() <= 1) {
        std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Spacegroup " << spacegroup << " does not contain at least one symmetry operation which is not the identity.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
      }

      topo.xray_asu().clear();
      for(int i = 0; i < spacegr.num_symops(); ++i) {
        int atom_pointer;
        _lineStream >> atom_pointer;
        if (_lineStream.fail()) {
          std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Cannot read ASU pointers " << (i+1) << ".";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }
        --atom_pointer;
        if (atom_pointer < 0 || atom_pointer >= int(topo.num_atoms())) {
          std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: ASU pointer " << (i+1) << " is out of range.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }
        DEBUG(9, "\t\tASU: " << atom_pointer);
        topo.xray_asu().push_back(atom_pointer);
      }
#else
      io::messages.add("In_Xrayresspec", "Compile with clipper support.", io::message::error);
#endif

      topo.xray_sym_restraints().clear();
      std::vector<std::string>::const_iterator it = buffer.begin() + 4,
              to = buffer.end() - 1;
      for(unsigned int line_nr = 6; it != to; ++it, ++line_nr) {
        std::string line(*it);
        if (line.length() < 17) {
         std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Line " << line_nr << " is too short.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }

        // the first 17 chars are ignored
        line.erase(line.begin(), line.begin() + 17);

        _lineStream.clear();
        _lineStream.str(line);

        int atom = 0;
        _lineStream >> atom;

        DEBUG(11, "\t" << atom);

        if (_lineStream.fail()) {
          std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Line " << line_nr << ": Cannot read atom.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }

        --atom;
        if (atom < 0 || atom >= int(topo.num_atoms())) {
          std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Line " << line_nr << ": Atom out of range.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }
        if (atom < int(topo.xray_asu()[0])) {
          std::ostringstream msg;
          msg << "XRAYSYMRESSPEC block: Line " << line_nr << ": Atom not in first ASU.";
          io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
          return;
        }
        const unsigned int atom_p = atom - topo.xray_asu()[0];
        for(unsigned int i = 1; i < topo.xray_asu().size(); ++i) {
          const unsigned int atom_img = topo.xray_asu()[i] + atom_p;
          if (atom_img >= topo.num_atoms()) {
            std::ostringstream msg;
            msg << "XRAYSYMRESSPEC block: Line " << line_nr << ": The image nr. "
                << i << " of atom " << (atom+1) << " is out of range.";
            io::messages.add("In_Xrayresspec", msg.str(), io::message::error);
            return;
          }
        }
        topo.xray_sym_restraints().push_back(atom);
      } // for atoms
    } // if block present
  } // XRAYSYMRESSPEC

  { // XRAYBFACTOROPTIMISATION

    buffer = m_block["XRAYBFACTOROPTIMISATION"];
    DEBUG(10, "XRAYBFACTOROPTIMISATION block : " << buffer.size());

    if (buffer.size()) {
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

      int step = 0, iter = 0;
      double grad = 0.0, bmin = 0.0, bmax = 0.0;
      // BFOPTS  BFOPTTI BFOPTTG BFOPTMN BFOPTMX
      _lineStream >> step >> iter >> grad >> bmin >> bmax;

      if (_lineStream.fail()) {
        io::messages.add("bad line in XRAYBFACTOROPTIMISATION block",
                "In_Xrayresspec", io::message::error);
        return;
      }
      if (step < 0) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTS has to be >= 0",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.bfactor.step = step;
      if (iter < 1) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTTI has to be >= 1",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.bfactor.terminate_iterations = iter;
      if (grad <= 0.0) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTTG has to be > 0.0",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.bfactor.terminate_gradient = grad;

      if (bmin < 0.0) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTMN has to be >= 0.0",
                "In_Xrayresspec", io::message::error);
        return;
      }
      if (bmax <= 0.0) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTMX has to be > 0.0",
                "In_Xrayresspec", io::message::error);
        return;
      }
      if (bmin >= bmax) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTMX has to be > BFOPTMN",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.bfactor.min = bmin;
      sim.param().xrayrest.bfactor.max = bmax;

      int num_group = 0;
      _lineStream >> num_group;
      if (_lineStream.fail() || num_group < 0) {
        io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTNG has to be >= 0",
                "In_Xrayresspec", io::message::error);
        return;
      }

      if (num_group == 0) {
        sim.param().xrayrest.bfactor.groups.resize(topo.num_atoms());
        for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
          sim.param().xrayrest.bfactor.groups[i].insert(i);
        }
      } else {
        sim.param().xrayrest.bfactor.groups.resize(num_group);
        for(unsigned int i = 0; i < unsigned(num_group); ++i) {
          int size = 0;
          _lineStream >> size;
          if (_lineStream.fail() || size <= 0) {
            io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTGS has to be > 0",
                    "In_Xrayresspec", io::message::error);
            return;
          }
          for (int j = 0; j < size; ++j) {
            int atom = 0;
            _lineStream >> atom;
            atom--;
            if (_lineStream.fail() || atom < 0 || atom >= int(topo.num_atoms())) {
              io::messages.add("XRAYBFACTOROPTIMISATION block: BFOPTGM, atom out of range.",
                      "In_Xrayresspec", io::message::error);
              return;
            }
            sim.param().xrayrest.bfactor.groups[i].insert(atom);
          }
        }
        // check for overlap
        for(unsigned int i = 0; i < unsigned(num_group); ++i) {
          std::set<unsigned int>::const_iterator it = sim.param().xrayrest.bfactor.groups[i].begin(),
                  to = sim.param().xrayrest.bfactor.groups[i].end();
          for(; it != to; ++it) {
            for(unsigned int j = i+1; j < unsigned(num_group); ++j) {
              if (sim.param().xrayrest.bfactor.groups[j].find(*it) !=
                  sim.param().xrayrest.bfactor.groups[j].end()) {
                std::ostringstream os;
                os << "XRAYBFACTOROPTIMISATION block: BFOPTGM, groups show overlap. "
                        << "Atom " << (*it)+1 << " is contained in group " << i+1
                        << " and " << j+1 << ".";
                io::messages.add(os.str(), "In_Xrayresspec", io::message::error);
              }
            } // for groups
          } // for members
        } // for groups
      } // if groups
    } // if buffer size
  } // XRAYBFACTOROPTIMISATION
  { // XRAYREPLICAEXCHANGE
    buffer = m_block["XRAYREPLICAEXCHANGE"];
    DEBUG(10, "XRAYREPLICAEXCHANGE block : " << buffer.size());
    if(buffer.size()){
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1,buffer.end() - 1,s));

      int interruptor = 0; //NTXRRE
      double min_value = 0.0, max_value = 0.0; //CXREEMN CXREEMX
      int energy_interruptor = 0;
      _lineStream >> interruptor >> min_value >> max_value >> energy_interruptor;

      if (_lineStream.fail()) {
        io::messages.add("bad line in XRAYREPLICAEXCHANGE block",
                "In_Xrayresspec", io::message::error);
        return;
      }
      switch (interruptor){
        case 0: sim.param().xrayrest.replica_exchange_parameters.switcher = simulation::replica_exchange_off;
        break;
        case 1: sim.param().xrayrest.replica_exchange_parameters.switcher = simulation::replica_exchange_force;
        break;
        case 2: sim.param().xrayrest.replica_exchange_parameters.switcher = simulation::replica_exchange_resolution;
        break;
        //std::cout << "interruptor" << interruptor<< std::endl;
        //std::cout << "switcher" << sim.param().xrayrest.replica_exchange_parameters.switcher << std::endl;
        default: io::messages.add("Forbidden value for NTXRRE",
                "In_Xrayresspec", io::message::error);
                return;
      }
      switch (energy_interruptor) {
        case 0: sim.param().xrayrest.replica_exchange_parameters.energy_switcher = simulation::energy_tot;
        break;
        case 1: sim.param().xrayrest.replica_exchange_parameters.energy_switcher = simulation::energy_phys;
        break;
        case 2: sim.param().xrayrest.replica_exchange_parameters.energy_switcher = simulation::energy_special;
        break;
        default: io::messages.add("Forbidden value for the energy interruptor",
                "In_Xrayresspec", io::message::error);
                return;
      }
      if (min_value >= max_value || min_value < 0.0){
        io::messages.add("Minimal and/or maximal force have/has absurd value(s)",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_min = min_value;
      sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_max = max_value;
    }
  } //XRAYREPLICAEXCHANGE
  { // XRAYOVERALLBFACTOR
    buffer = m_block["XRAYOVERALLBFACTOR"];
    DEBUG(10, "XRAYOVERALLBFACTOR block : " << buffer.size());
    if(buffer.size()){
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1,buffer.end() - 1,s));

      double overall_Bfac = 0.0;
      int overall_Bfac_interruptor = 0;
      _lineStream >> overall_Bfac >> overall_Bfac_interruptor;

      if (_lineStream.fail()) {
        io::messages.add("bad line in XRAYOVERALLBFACTOR block",
                "In_Xrayresspec", io::message::error);
        return;
      }
      switch (overall_Bfac_interruptor){
        case 0: sim.param().xrayrest.overall_bfactor.B_overall_switcher = simulation::B_overall_off;
        break;
        case 1: sim.param().xrayrest.overall_bfactor.B_overall_switcher = simulation::B_overall_on;
        break;
        default: io::messages.add("Forbidden value for XROBF",
                "In_Xrayresspec", io::message::error);
                return;
      }
      sim.param().xrayrest.overall_bfactor.init = overall_Bfac;
    }
  } //XRAYOVERALLBFACTOR
  { //XRAYSFCALC
    buffer = m_block["XRAYSFCALC"];
    if(buffer.size()){
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1,buffer.end() - 1,s));

      double sf_tolerance = 0.0; //SFCTOL
      unsigned int sf_constant = 0; //SFCST
      _lineStream >> sf_tolerance >> sf_constant;

      if (_lineStream.fail()) {
        io::messages.add("bad line in XRAYSFCALC block",
                "In_Xrayresspec", io::message::error);
        return;
      }

      if (sf_tolerance < 0.0){
        io::messages.add("XRAYSFCALC block: tolerance has negative value",
                "In_Xrayresspec", io::message::error);
        return;
      }
      sim.param().xrayrest.structure_factor_calculation.atom_move_tolerance = sf_tolerance;
      sim.param().xrayrest.structure_factor_calculation.steps_nb_constant = sf_constant;
    }
  }//XRAYSFCALC
  if (!quiet) {
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
                os << "\tx-ray restraints OFF\n";
        break;
      case simulation::xrayrest_inst :
                os << "\tx-ray instantaneous restraints ON\n";
        break;
      case simulation::xrayrest_avg :
                os << "\tx-ray time-averaged restraints ON\n";
        break;
      case simulation::xrayrest_biq :
                os << "\tx-ray biquadratic instantaneous/time-averaged restraints ON\n";
        break;
    }

    switch (sim.param().xrayrest.symrest) {
      case simulation::xray_symrest_off:
        os << "\tsymmetry restraints OFF\n";
        break;
      case simulation::xray_symrest_ind:
        os << "\tsymmetry restraints on individual atom positions\n";
        break;
      case simulation::xray_symrest_constr:
        os << "\tsymmetry constraints on individual atom positions\n";
        break;
    }
    os << "END\n";
  }
}


