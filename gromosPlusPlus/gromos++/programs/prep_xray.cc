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
 * @file prep_xray.cc
 * Creates input-file for xray restraining
 */

/**
 * @page programs Program Documentation
 *
 * @anchor prep_xray
 * @section prep_xray Creates input-file for crystallographic restraining
 * @author @ref ff @ref ns
 * @date 3-2-2009
 *
 * Program prep_xray converts a crystallographic information file (CIF) containing
 * reflection data into a GROMOS X-ray restraints specification file. Using a
 * mapping file (\@map) it writes out the element names of the solute and solvent
 * atoms according to their integer atom codes. The atoms' B-factors and occupancies are read from a
 * special file (\@bfactor) if requested or defaulted to @f$ 0.01 \mathrm{nm}^2 @f$ and 100%.
 * The reflection list can be filtered according to the given resolution range. 
 * If no resolution is given, it is determined automatically.
 * A random set of amplitudes is created for the computation of the free R factor.
 * 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@cif</td><td>&lt;crystallographic information file&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt; @ref gio::InIACElementNameMapping "file with IAC-to-elementname mapping" &gt; </td></tr>
 * <tr><td> \@spacegroup</td><td>&lt;spacegroup in Hermann-Maauguin format&gt; </td></tr>
 * <tr><td> \@cell</td><td>&lt;cell in form: a b c alpha beta gamma&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution, from and to&gt; </td></tr>
 * <tr><td> \[\@filter</td><td>&lt;filter off small structure factor amplitudes&gt;]</td></tr>
 * <tr><td> \@bfactor</td><td>&lt;@ref gio::InBFactorOccupancy "a B factor and occupancies file"&gt;</td></tr>
 * <tr><td> \@symmetrise</td><td>&lt;apply symmetry operations to relection list &gt;</td></tr>
 * <tr><td> \@rfree</td><td>&lt;percentage of reflections used for r_free calculation, random number seed&gt;</td></tr>
 * <tr><td>[\@factor</td><td>&lt;convert length unit to Angstrom&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  prep_xray
    @topo          ex.top
    @cif           ex.cif
    @map           ex.map
    @spacegroup    P 21 21 21
    @cell          5.00 5.10 5.20 90.0 90.0 90.0
    @resolution    0.3 0.15
    @filter        2.0
    @bfactor       ex.boc
    @rfree         10 1234
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <limits>
#include <ios>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InCIF.h"
#include "../src/gio/InIACElementNameMapping.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/gmath/Stat.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace std;

// Additional Clipper Headers
#include "../config.h"
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>

int main(int argc, char *argv[]) {
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@cif            <cristallographic information file>\n";
  usage += "\t@map            <gromos atom code to element mapping file>\n";
  usage += "\t@spacegroup     <spacegroup in hermann-mauguin form>\n";
  usage += "\t@cell           <cell in form: a b c alpha beta gamma>\n";
  usage += "\t@resolution     <scattering resolution min max>\n";
  usage += "\t@bfactor        <B-factor and occupancy file>\n";
  usage += "\t@symmetrise     <apply symmetry operation to reflections>\n";
  usage += "\t@rfree          <percentage of HKLs used for R-free, seed>\n";
  usage += "\t[@filter        <filter off small structure factor amplitudes>]\n";
  usage += "\t[@factor     <convert length unit to Angstrom. default: 10.0>]\n";


  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "cif" << "map" << "spacegroup" << "cell" << "resolution"
          << "bfactor" << "symmetrise" << "rfree" << "filter" << "factor";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    double factor = args.getValue<double>("factor", false, 10.0);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //READ IN MANUAL VALUES
    // Get Spacegroup Data
    string spgr;
    {
      Arguments::const_iterator iter = args.lower_bound("spacegroup");
      Arguments::const_iterator to = args.upper_bound("spacegroup");
      spgr = iter->second;
      for (++iter; iter != to; ++iter) {
        spgr += " ";
        spgr += iter->second;
      }
    }

    // initialize the spacegroup
    clipper::Spacegroup spacegroup;
    try {
      spacegroup.init(clipper::Spgr_descr(spgr, clipper::Spgr_descr::HM));
    } catch(clipper::Message_fatal & msg) {
      throw gromos::Exception(argv[0], "Invalid spacegroup: " + msg.text());
    }

    //Read in cell
    vector<float> cell;
    {
      Arguments::const_iterator iter = args.lower_bound("cell");
      Arguments::const_iterator to = args.upper_bound("cell");
      int i = 0;
      for (; iter != to; ++iter, ++i) {
        float tcell;
        std::istringstream is(iter->second);
        if (!(is >> tcell)) {
          throw gromos::Exception(argv[0],
                  "Cell parameters not numeric");
        }
        cell.push_back(tcell);
      }
      if (i != 6) {
        throw gromos::Exception(argv[0],
                "Not enough cell parameters");
      }
    }
    // create the clipper cell
    clipper::Cell clipperCell(
            clipper::Cell_descr(cell[0]*factor, cell[1]*factor, cell[2]*factor,
                cell[3], cell[4], cell[5]));

    //read in scattering resolution
    bool have_reso = false;
    if (args.count("resolution") == 2) {
      have_reso = true;
    } else {
      throw gromos::Exception(argv[0],
              "Please give a resolution range or let the program determine it");
    }
    vector<double> resoarg = args.getValues<double>("resolution", 2, false,
          Arguments::Default<double>() << 0.0 << numeric_limits<double>::max());
    double reso_min = resoarg[0];
    double reso_max = resoarg[1];

    bool has_filter = false;
    if (args.count("filter") == 1) {
      has_filter = true;
    }
    double filter = args.getValue<double>("filter", false, 0.0);

    // swap them if someone is confused by the fact that the low numeric
    // value is actually a high resolution
    if (have_reso && reso_min <= reso_max) {
      const double tmp = reso_min;
      reso_min = reso_max;
      reso_max = tmp;
    }

    // Read in cif-file for generating the structure factor list (using clipper)
    InCIF ciffile(args["cif"]);
    vector<CIFData> cifdata = ciffile.getData();

    InBFactorOccupancy bocfile(args["bfactor"]);
    vector<BFactorOccupancyData> bocdata = bocfile.getData();

    // Generate Mapping
    InIACElementNameMapping mapfile(args["map"]);
    map<int, string> gacmapping = mapfile.getData();
    vector<string> element;
    vector<BFactorOccupancyData> boc;
    unsigned int c = 0;
    for (int j = 0; j < sys.numMolecules(); ++j) {
      for (int i = 0; i < sys.mol(j).numAtoms(); ++i, ++c) {
        element.push_back(gacmapping[sys.mol(j).topology().atom(i).iac()]);
        if (c >= bocdata.size())
          throw gromos::Exception(argv[0], "Not enough B factors and occupancies given.");
        boc.push_back(bocdata[c]);
      }
    }
    // Generate Solvent Mapping
    vector<string> solvent_element;
    vector<BFactorOccupancyData> solvent_boc;
    for (int j=0; j<sys.numSolvents(); ++j){
      for (int i = 0; i < sys.sol(j).topology().numAtoms(); ++i,++c) {
        solvent_element.push_back(gacmapping[sys.sol(j).topology().atom(i).iac()]);
        if (c >= bocdata.size())
          throw gromos::Exception(argv[0], "Not enough B factors and occupancies for solvent given. Add them after the solute.");
        solvent_boc.push_back(bocdata[c]);
      }
    }

    if (has_filter) {
      gmath::Stat<double> sfstat;
      for (unsigned int i = 0; i < cifdata.size(); ++i)
        sfstat.addval(cifdata[i].f_obs);

      const double threshold = sfstat.ave() - filter * sfstat.rmsd();
      cerr << "Threshold for structure factor: " << threshold << endl;
      
      vector<CIFData> filtered;
      for(unsigned int i = 0; i < cifdata.size(); ++i) {
        if (cifdata[i].f_obs >= threshold) {
          filtered.push_back(cifdata[i]);
        }
      }
      cerr << (cifdata.size() - filtered.size()) << " structure factor amplitudes have been filtered away." << endl;
      cifdata = filtered;
    }

    if (args.count("symmetrise") >= 0) {;
      const unsigned int orig_size = cifdata.size();
      for (unsigned int i = 0; i < orig_size; ++i) {
        const clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
        // loop over symmetry operations
        for(int j = 1; j < spacegroup.num_symops(); ++j) {
          const clipper::HKL & hkl_img = hkl.transform(spacegroup.symop(j));
          if (hkl_img == hkl)
            continue;

          // store new reflection
          CIFData img;
          img.h = hkl_img.h(); img.k = hkl_img.k(); img.l = hkl_img.l();
          img.f_obs = cifdata[i].f_obs;
          img.stddev_f_obs = cifdata[i].stddev_f_obs;
          cifdata.push_back(img);
        } // loop over symmetry operations
      } // loop over reflections
    } // symmetrise

    vector<double> rfreearg = args.getValues<double>("rfree", 2, false,
          Arguments::Default<double>() << 0.0 << 1234);
    double r_free = rfreearg[0];
    int seed = int(rfreearg[1]);
    if (r_free < 0.0 || r_free > 100.0) {
      throw gromos::Exception(argv[0], "@rfree: percentage must be 0..100");
    }

    vector<CIFData> restrCif;
    vector<CIFData> freeCif;
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    for(unsigned int i = 0; i < cifdata.size(); ++i) {
      // flip coin
      double coin = gsl_rng_uniform(rng) * 100.0;
      if (coin < r_free)
        freeCif.push_back(cifdata[i]);
      else
        restrCif.push_back(cifdata[i]);
    }

    

    // OUTPUT
    // DUMMY TITLE
    cout << "TITLE\n";
    cout << "Xray specification file\n";
    cout << "END\n";

    // XRAYRESSPEC
    cout << "XRAYRESSPEC\n";
    cout << "#" << setw(4) << "H" << setw(5) << "K" << setw(5) << "L";
    cout << setw(10) << "SF" << setw(12) << "STDDEV_SF\n";

    int filtered = 0;
    for (unsigned int i = 0; i < restrCif.size(); ++i) {
      clipper::HKL hkl(restrCif[i].h, restrCif[i].k, restrCif[i].l);
      double hkl_resolution = sqrt(1.0 / hkl.invresolsq(clipperCell)) / factor;
      if (have_reso) {
        if (hkl_resolution < reso_max || hkl_resolution > reso_min) {
          ++filtered;
          continue;
        }
      } else {
        // determine maximal resolution (minimal numerical number)
        reso_max = std::min(reso_max, hkl_resolution);
      }
      cout << setw(5) << restrCif[i].h;
      cout << setw(5) << restrCif[i].k;
      cout << setw(5) << restrCif[i].l;
      cout << setw(10) << restrCif[i].f_obs;
      cout << setw(11) << restrCif[i].stddev_f_obs << "\n";
    }
    cout << "END\n";

    // XRAYRFREESPEC
    cout << "XRAYRFREESPEC\n";
    cout << "#" << setw(4) << "H" << setw(5) << "K" << setw(5) << "L";
    cout << setw(10) << "SF" << setw(12) << "STDDEV_SF\n";

    for (unsigned int i = 0; i < freeCif.size(); ++i) {
      clipper::HKL hkl(freeCif[i].h, freeCif[i].k, freeCif[i].l);
      double hkl_resolution = sqrt(1.0 / hkl.invresolsq(clipperCell)) / factor;
      if (have_reso) {
        if (hkl_resolution < reso_max || hkl_resolution > reso_min) {
          ++filtered;
          continue;
        }
      } else {
        // determine maximal resolution (minimal numerical number)
        reso_max = std::min(reso_max, hkl_resolution);
      }
      cout << setw(5) << freeCif[i].h;
      cout << setw(5) << freeCif[i].k;
      cout << setw(5) << freeCif[i].l;
      cout << setw(10) << freeCif[i].f_obs;
      cout << setw(11) << freeCif[i].stddev_f_obs << "\n";
    }
    cout << "END\n";

    if (filtered)
      cerr << filtered << " reflections have been filtered away." << endl;

    //XRAYELEMENTSPEC
    cout << "XRAYELEMENTSPEC\n";
    cout << "# ELEMENT[1..N]\n";
    for (unsigned int i = 0; i < element.size(); ++i) {
      cout << setw(2) << element[i] << " ";
      if ((i + 1) % 20 == 0) {
        cout << "\n";
      }
    }
    cout << "\nEND\n";

    //XRAYSOLVELEMENTSPEC
    cout << "XRAYSOLVELEMENTSPEC\n";
    cout << "# ELEMENT[1..N]\n";
    for (unsigned int i = 0; i < solvent_element.size(); ++i) {
      cout << setw(2) << solvent_element[i] << " ";
      if ((i + 1) % 20 == 0) {
        cout << "\n";
      }
    }
    cout << "\nEND\n";

    //XRAYBFOCCSPEC
    cout.precision(6);
    cout << "XRAYBFOCCSPEC\n";
    cout << "# BF    OCC\n";
    for(unsigned int i = 0; i < boc.size(); ++i) {
      cout << setw(15) << boc[i].b_factor << setw(10) << boc[i].occupancy << "\n";
    }
    cout << "END\n";
    //XRAYSOLVBFOCCSPEC
    cout << "XRAYSOLVBFOCCSPEC\n";
    cout << "# BF    OCC\n";
    for(unsigned int i = 0; i < solvent_boc.size(); ++i) {
      cout << setw(15) << solvent_boc[i].b_factor << setw(10) << solvent_boc[i].occupancy << "\n";
    }
    cout << "END\n";
    

    // XRAYRESPARA
    cout.precision(4);
    cout << "XRAYRESPARA\n";
    cout << "#" << setw(15) << "SPGR\n";
    cout << setw(15) << spgr << "\n";
    cout << "#" << setw(9) << "RESO" << setw(10) << "TOANG" << "\n";
    cout << setw(10) << reso_max << setw(10) << factor<< "\n";
    cout << "END\n";

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
