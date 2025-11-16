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
 * @file r_factor.cc
 * calculates r factors
 */

/**
 * @page programs Program Documentation
 *
 * @anchor r_real_factor
 * @section r_real_factor calculates the crystallographic real-space residual (real-space R factors)
 * @author @ref ns
 * @date 14.02.2010
 *
 * Program r_real_factor calculates two electron densities. One (@f$\rho_\mathrm{calc}@f$) from the atomic positions
 * and a second (@f$\rho_\mathrm{obs}@f$) from the structure factor amplitudes and calculated phases.
 * Only the atoms given by the @ref AtomSpecifier \@atomssf are considered for 
 * the structure factor calculation.
 *
 * The real space residual
 * @f[ R = \frac{\sum\alpha\rho_\mathrm{obs} + \beta - \rho_\mathrm{calc}}{\sum\alpha\rho_\mathrm{obs} + \beta + \rho_\mathrm{calc}} @f]
 * is calculated for every residue. Summation is only carried out over the extent
 * of the atoms contained in the @ref AtomSpecifier \@atomsr
 *
 * The atoms' IAC are mapped to their element names according to the rules given
 * in the \@map file. The atoms' B-factors and occupancies are read from a
 * special file (\@bfactor) if requested or defaulted to @f$ 0.01 \mathrm{nm}^2 @f$ and 100%.
 * The electron densities are calculated to the given resolution (\@resultion) while
 * the cell information is calculated from the system's box.
 * Symmetry operations are taken into account by specifing a (\@spacegroup).
 * Make sure you only give asymetric unit when using \@spacegroup.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomssf</td><td>&lt;@ref AtomSpecifier atoms to consider for structure_factor&gt; </td></tr>
 * <tr><td> \@atomsr</td><td>&lt;@ref AtomSpecifier atoms to consider for R value&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;file with IAC-to-elementname mapping&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;file with experimental B-factors and occupancies&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution&gt; </td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup in Hermann-Mauguin format, default: P 1&gt;]</td></tr>
 * <tr><td> \@cif</td><td>&lt;cristallographic information file&gt; </td></tr>
 * <tr><td>[\@factor</td><td>&lt;convert length unit to Angstrom&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 r_real_factor
    @topo       ex.top
    @time       0 0.1
    @atomssf    a:a
    @atomsr     1:C,O,N,CA
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
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Physics.h"
#include "../src/gio/InIACElementNameMapping.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/gio/InCIF.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/debug.h"

// Additional Clipper Headers
#include "../config.h"
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <set>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace std;
using namespace gmath;
using namespace bound;

/**
 * fit two electron densities on top of each other. This is done using a
 * linear regression such that
 * @f[ \rho_1 = \alpha + \beta\rho_2 @f]
 *
 * @param[in] rho1 the first electron density
 * @param[in] rho2 the second electron density
 * @param[in] points the set of grid point to consider
 * @param[out] slope slope @f$\beta@f$ of the linear regression
 * @param[out] intercept intercept @f$\alpha@f$  of the linrar regression
 */
void fit_rho(
        const clipper::Xmap<clipper::ftype64> & rho1,
        const clipper::Xmap<clipper::ftype64> & rho2,
        std::set<int> & points,
        double & slope, double & intercept) {
  DEBUG(6, "fitting electron densities");
  // if the set is empty just fill it with all grid points
  if (points.empty()) {
    DEBUG(10, "\tfitting whole grid.");
    for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = rho1.first();
          !ix.last(); ix.next()) {
      points.insert(ix.index());
    }
  }

  double sum_xy = 0.0, sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0;
  const double n = points.size();
  for(std::set<int>::const_iterator it = points.begin(), to = points.end(); it != to; ++it) {
    // get the data
    const double & x = rho2.get_data(*it);
    const double & y = rho1.get_data(*it);
    // calculate the suns
    sum_xy += x*y;
    sum_x += x;
    sum_y += y;
    sum_xx += x*x;
  }
  const double mean_x = sum_x / n;
  DEBUG(9, "mean rho2: " << mean_x);
  const double mean_y = sum_y / n;
  DEBUG(9, "mean rho1: " << mean_y);
  slope = (sum_xy - sum_x*mean_y) / (sum_xx - sum_x*mean_x);
  intercept = mean_y - slope * mean_x;
  DEBUG(9, "slope = " << slope << " intercept = " << intercept);
}


int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "map" << "atomssf" << "time" << "bfactor"
          << "resolution" << "spacegroup" << "cif" << "atomsr" << "factor";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atomssf     <atomspecifier: atoms to consider for structure_factor>\n";
  usage += "\t@atomsr      <atomspecifier: atom used for R real calculation>\n";
  usage += "\t@traj        <trajectory files>\n";
  usage += "\t@map         <IAC-to-ElementName map-file>\n";
  usage += "\t[@bfactor    <experimental B-factors>]\n";
  usage += "\t@resolution  <scattering resolution>\n";
  usage += "\t[@spacegroup <spacegroup in Hermann-Maugin format, default: P 1>]\n";
  usage += "\t@cif         <cristallographic information file>\n";
  usage += "\t[@factor     <convert length unit to Angstrom. default: 10.0>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);



  try {
    Arguments args(argc, argv, knowns, usage);

    double factor = args.getValue<double>("factor", false, 10.0);

    // Hardcoded B-factor conversion factor.
    const double sqpi2=(M_PI*M_PI*8.0);

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

    AtomSpecifier atomsr(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsr");
      Arguments::const_iterator to = args.upper_bound("atomsr");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atomsr.addSpecifier(spec);
      }
    }
    if (atomsr.size() == 0)
      throw gromos::Exception(argv[0], "No R-factor-atoms specified!");

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

    // Read in cif-file for generating the structure factor list (using clipper)
    InCIF ciffile(args["cif"]);
    vector<CIFData> cifdata = ciffile.getData();

    //===========================
    // loop over all trajectories
    InG96 ic;

    cout << "#" << setw(7) << "residue" << setw(15) << "R real" << endl;
    unsigned int frame_number = 0;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        ++frame_number;
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
        for(int i = 0; i < atomsr.size(); ++i) {
          atomsr.pos(i) = atomsr.pos(i) - pbc->nearestImage(atomsr.pos(i), centre, sys.box()) + centre;
        }

        // create the cell
        const double a = sys.box().K().abs();
        const double b = sys.box().L().abs();
        const double c = sys.box().M().abs();
        const double alpha = sys.box().alpha();
        const double beta = sys.box().beta();
        const double gamma = sys.box().gamma();

        if (!a || !b || !c || !alpha || !beta || !gamma)
          throw gromos::Exception(argv[0], "Box has zero volume!");

        clipper::Cell_descr cellinit(a * factor, b * factor, c * factor,
                alpha, beta, gamma);
        clipper::CCell cell(spgr, clipper::String("base cell"), clipper::Cell(cellinit));

        // create the resolutions and corresponding lattice
        clipper::CResolution reso(cell, clipper::String("base reso"), clipper::Resolution(resolution * factor));
        clipper::CHKL_info hkls(reso, clipper::String("base hkls"), true);
        clipper::CHKL_data<clipper::data64::F_phi> fphi(hkls);
        clipper::CHKL_data<clipper::data64::F_phi> fphi_print(hkls);

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
            atm.set_u_iso(bfoc[atom_index].b_factor * factor*factor / sqpi2);
          } else {
            atm.set_occupancy(1.0);
            atm.set_u_iso(1.0 / sqpi2);
          }
          atm.set_element(gacmapping[calcatoms.iac(i)]);
          if (atm.occupancy() != 0.0) {
            atomvec.push_back(atm);
          }
        }
        clipper::Atom_list atoms(atomvec);

        // Calculate structure factors
        clipper::SFcalc_iso_fft<double> sfc;
        sfc(fphi, atoms);

        
        fphi_print = std::complex<double>(0.0,0.0);
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          const clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
          fphi_print.set_data(hkl, clipper::data64::F_phi(cifdata[i].f_obs, fphi[hkl].phi()));
        }

        // calculate densities
        const clipper::Grid_sampling grid_obs(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
        clipper::Xmap<clipper::ftype64> rho_obs(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid_obs);
        rho_obs.fft_from(fphi_print);
        const clipper::Grid_sampling grid_calc(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
        clipper::Xmap<clipper::ftype64> rho_calc(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid_calc);
        rho_calc.fft_from(fphi);

        // loop over residues
        map<int,set<int> > residues;
        for (int i = 0; i < atomsr.size(); i++) {
          int resnum = atomsr.resnum(i);
          for(int m = 0; m < atomsr.mol(i); ++m)
            resnum += sys.mol(m).topology().numRes();
          residues[resnum].insert(i);
        }
        
        map<int,set<int> >::const_iterator res_it = residues.begin(),
                res_to = residues.end();
        const double cutoff = 0.25 * factor;
        const double cutoff2 = cutoff * cutoff;
        clipper::Grid_range gd(cell, grid_obs, cutoff);
        for(; res_it != res_to; ++res_it) {
          set<int> points;
          set<int>::const_iterator atom_it = res_it->second.begin(),
                  atom_to = res_it->second.end();
          for (; atom_it != atom_to; ++atom_it) {
            gmath::Vec r = atomsr.pos(*atom_it);
            r *= factor;
            const clipper::Coord_orth atom_pos(r[0], r[1], r[2]);

            // determine grad range of atom
            const clipper::Coord_frac uvw = atom_pos.coord_frac(cell);
            const clipper::Coord_grid g0 = uvw.coord_grid(grid_obs) + gd.min();
            const clipper::Coord_grid g1 = uvw.coord_grid(grid_obs) + gd.max();

            // loop over atom's grid
            clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
            i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_obs, g0);
            for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
              for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
                for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                  // check if this cell is within the cutoff
                  if ((iw.coord_orth() - atom_pos).lengthsq() < cutoff2) {
                    points.insert(iw.index());
                  } // if within cutoff
                }
              }
            } // loop over grid
          } // loop over atoms
          // the fitting parameters
          double a, b;
          fit_rho(rho_calc, rho_obs, points, a, b);
          double sum_obs_m_calc = 0.0;
          double sum_obs_p_calc = 0.0;
          // loop over grid points
          for (std::set<int>::const_iterator it = points.begin(), to = points.end();
                  it != to; ++it) {
            // bring obs on the same scale as calc
            const double obs = a * rho_obs.get_data(*it) + b;
            const double calc = rho_calc.get_data(*it);
            // for R-real
            sum_obs_m_calc += fabs(obs - calc);
            sum_obs_p_calc += fabs(obs + calc);
          }
          const double r_real = sum_obs_m_calc / sum_obs_p_calc;

          cout.precision(8);
          cout << setw(8) << (res_it->first + 1) << setw(15) << r_real << endl;
        } // loop over residue
        cout << endl;
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
