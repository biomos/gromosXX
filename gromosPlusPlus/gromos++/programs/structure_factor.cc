/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file structure_factor.cc
 * calculates structure factors
 */
/**
 * @page programs Program Documentation
 *
 * @anchor structure_factor
 * @section structure_factor calculates structure factors
 * @author @ref ns ff
 * @date 8.4.2009
 *
 * Program structure_factor calculates crystallographic structure-factor amplitudes
 * and phases from a given trajectory. Only the atoms given by the @ref AtomSpecifier
 * \@atomssf are considered for the calculation. The atoms' IAC are mapped to their
 * element names according to the rules given in the \@map file. The atoms' B-factors
 * and occupancies are read from a special file (\@bfactor) if requested or defaulted
 * to @f$ 0.01 \mathrm{nm}^2 @f$ and 100%.
 * Structure-factor amplitudes are calculated to the given resolution (\@resultion) while
 * the cell information is calculated from the system's box.
 * Symmetry operations are taken into account by specifying a space group (\@spacegroup).
 * When using \@spacegroup, make sure only the asymmetric unit is given.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomssf</td><td>&lt;@ref AtomSpecifier atoms to consider for structure_factor&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;@ref gio::InIACElementNameMapping "file with IAC-to-elementname mapping" &gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;@ref gio::InBFactorOccupancy "file with experimental B-factors and occupancies"&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution [nm]&gt; </td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup in Hermann-Mauguin format, default: P 1&gt;]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;convert length unit to Angstrom&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  structure_factor
    @topo       ex.top
    @time       0 0.1
    @atomssf    1:CA
    @traj       ex.tr
    @map        ex.map
    @bfactor    ex.bfc
    @resolution 0.1
    @spacegroup P 21 21 21
    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Physics.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InIACElementNameMapping.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/utils/debug.h"

// Additional Clipper Headers
#include "../config.h"
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace std;
using namespace gmath;
using namespace bound;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "map" << "atomssf" << "time" << "bfactor"
          << "resolution" << "spacegroup" << "factor";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atomssf     <atomspecifier: atoms to consider for structure_factor>\n";
  usage += "\t@traj        <trajectory files>\n";
  usage += "\t@map         <IAC-to-ElementName map-file>\n";
  usage += "\t[@bfactor    <experimental B-factors>]\n";
  usage += "\t@resolution  <scattering resolution>\n";
  usage += "\t[@spacegroup <spacegroup in Hermann-Maugin format, default: P 1>]\n";
  usage += "\t[@factor     <convert length unit to Angstrom. default: 10.0>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);

  try {
    Arguments args(argc, argv, knowns, usage);

    // Hardcoded B-factor conversion factor.
    const double sqpi2=(M_PI*M_PI*8.0);

    double factor = args.getValue<double>("factor", false, 10.0);

    // Get Spacegroup Data or default to no symmetry (P 1)
    std::string spgrdata("P 1");
    {
      Arguments::const_iterator iter = args.lower_bound("spacegroup");
      Arguments::const_iterator to = args.upper_bound("spacegroup");
      if (iter != to) {
        spgrdata = iter->second;
        for (++iter; iter != to; ++iter) {
          spgrdata += " ";
          spgrdata += iter->second;
        }
      }
    }
    // initialize the spacegroup
    std::auto_ptr<clipper::Spgr_descr> spgrinit;
    try {
     spgrinit = std::auto_ptr<clipper::Spgr_descr>(new clipper::Spgr_descr(spgrdata, clipper::Spgr_descr::HM));
    } catch(clipper::Message_fatal & msg) {
      throw gromos::Exception(argv[0], "Invalid spacegroup: " + msg.text());
    }
    clipper::CSpacegroup spgr(clipper::String("base spgr"), clipper::Spacegroup(*spgrinit));

    // Get resolution as a double
    double resolution = args.getValue<double>("resolution", true);

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    // System for calculation
    System sys(it.system());
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    AtomSpecifier calcatoms(sys);
    //get structure_factor atoms
    {
      Arguments::const_iterator iter = args.lower_bound("atomssf");
      Arguments::const_iterator to = args.upper_bound("atomssf");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        calcatoms.addSpecifier(spec);
      }
    }
    if (calcatoms.size() == 0)
      throw gromos::Exception(argv[0], "No structure_factor-atoms specified!");
    
    //Get gac-to-ele mapping
    InIACElementNameMapping mapfile(args["map"]);
    std::map<int, std::string> gacmapping = mapfile.getData();

    // Get experimental Bfactors and occupancy
    std::vector<BFactorOccupancyData> bfoc;
    bool has_bfactor = false;
    if (args.count("bfactor") == 1) {
      InBFactorOccupancy bfac_file(args["bfactor"]);
      bfoc = bfac_file.getData();
      has_bfactor = true;
    }

    //===========================
    // loop over all trajectories
    InG96 ic;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        if (!sys.hasPos)
          throw gromos::Exception(argv[0],
                "Unable to read POSITION(RED) block from "
                "trajectory file.");
        
        if (!sys.hasBox)
          throw gromos::Exception(argv[0],
                "Cannot calculate structure factors without a box.");

        // get the centre of the box
        const Vec centre = (sys.box().K() * 0.5) + (sys.box().L() * 0.5) +
                (sys.box().M() * 0.5);

        // put atom into positive box
        for(int i = 0; i < calcatoms.size(); ++i) {
          calcatoms.pos(i) = calcatoms.pos(i) - pbc->nearestImage(calcatoms.pos(i), centre, sys.box()) + centre;
        }

        // create the cell
        clipper::Cell_descr cellinit(sys.box().K().abs() * factor, sys.box().L().abs() * factor, sys.box().M().abs() * factor,
                sys.box().alpha(), sys.box().beta(), sys.box().gamma());
        clipper::CCell cell(spgr, clipper::String("base cell"), clipper::Cell(cellinit));

        // create the resolutions and corresponding lattice
        clipper::CResolution reso(cell, clipper::String("base reso"), clipper::Resolution(resolution * factor));
        clipper::CHKL_info hkls(reso, clipper::String("base hkls"), true);
        clipper::CHKL_data<clipper::data64::F_phi> fphi(hkls);

        // Fill Clipper Atom list
        // we do this insight the loop due to solvent molecules!
        std::vector<clipper::Atom> atomvec;
        for (int i = 0; i < calcatoms.size(); i++) {
          clipper::Atom atm;
          // convert to angstrom
          atm.set_coord_orth(clipper::Coord_orth(
                  calcatoms.pos(i)[0] * factor,
                  calcatoms.pos(i)[1] * factor,
                  calcatoms.pos(i)[2] * factor));

          if (has_bfactor) {
            const unsigned int atom_index = calcatoms.gromosAtom(i);
            if (atom_index >= bfoc.size()) {
              throw gromos::Exception("structre_factor", "Not enough B-factors given");
            }
            atm.set_occupancy(bfoc[atom_index].occupancy);
            // convert to Angstrom^2
            atm.set_u_iso(bfoc[atom_index].b_factor * factor * factor / sqpi2);
          } else {
            atm.set_occupancy(1.0);
            atm.set_u_iso(1.0 / sqpi2);
          }
          atm.set_element(gacmapping[calcatoms.iac(i)]);
          atomvec.push_back(atm);
        }
        clipper::Atom_list atoms(atomvec);

        // Calculate structure factors
        clipper::SFcalc_iso_fft<double> sfc;
        sfc(fphi, atoms);

        cout << "# time: " << time << endl;
        cout << "loop_" << endl
                << "_refln.index_h" << endl
                << "_refln.index_k" << endl
                << "_refln.index_l" << endl
                << "_refln.F_meas_au" << endl
                << "_refln.phase_meas" << endl;

        for (clipper::HKL_info::HKL_reference_index ih = fphi.first_data(); !ih.last(); fphi.next_data(ih)) {
          cout << setw(6) << ih.hkl().h()
                  << setw(6) << ih.hkl().k()
                  << setw(6) << ih.hkl().l()
                  << setw(15) << fphi[ih].f()
                  << setw(15) << fphi[ih].phi() * 180.0 / M_PI << "\n";
        }

      } // while frames in file
    } // for traj
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
#else
int main(int argc, char **argv) {
  cerr << "You have to compile GROMOS++ with CCP4"
              " and clipper libraries in order to use this program.\n"
              "Use --with-ccp4 and --with-clipper for configuration."
          << endl;
  return 1;
}
#endif
