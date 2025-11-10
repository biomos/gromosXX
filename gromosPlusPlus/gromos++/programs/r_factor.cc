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
 * @anchor r_factor
 * @section r_factor calculates crystallographic R factors
 * @author @ref ns ff
 * @date 8.4.2009
 *
 * Program r_factor calculates crystallographic structure-factor amplitudes
 * and phases from a given trajectory and compares them to experimental values.
 * Only the atoms given by the @ref AtomSpecifier
 * \@atomssf are considered for the calculation. The atoms' IAC are mapped to their
 * element names according to the rules given in the \@map file. The atoms' B-factors
 * and occupancies are read from a special file (\@bfactor) if requested or defaulted
 * to @f$ 0.01 \mathrm{nm}^2 @f$ and 100%.
 * Structure factors are calculated to the given resolution (\@resultion) while
 * the cell information is calculated from the system's box.
 * Symmetry operations are taken into account by specifying a space group (\@spacegroup).
 * Make sure you only give asymmetric unit when using \@spacegroup.
 * The program can write the electron density (a @f$ 2|F_\mathrm{obs}| - |F_\mathrm{calc}| @f$ map)
 * to special files (FRAME_DENSITY_00001.ccp4), if requested (\@density flag).
 *
 * A bulk solvent correction can be applied if \@solvent is given. In a first step
 * a solvent mask is determined. Therefore the parameters @f$ r_\mathrm{vdW} @f$,
 * @f$ r_\mathrm{ion} @f$, @f$ r_\mathrm{shrink} @f$ and the IAC of the water oxygen
 * have to be provided. The occupied space is determined by the van-der-Waals
 * radius of the atoms plus a probe radius. The van-der-Waals radius of an atom
 * is calculated as half of the distance where the Lennard Jones potential energy
 * of the atom/water-oxygen interaction reaches its minimum. The probe radius is
 * either taken as @f$ r_\mathrm{vdW} @f$ (for neutral atoms) or @f$ r_\mathrm{ion} @f$
 * (for charged atoms). The occupied space is shrinked by @f$ r_\mathrm{shrink} @f$.
 * The structure factor is calculated as
 * @f[ F_\mathrm{calc} = F_\mathrm{model} + \rho\exp(-B\sin(\theta)^2/\lambda^2)\mathcal{F}(\mathbf{M}) \mathrm{.} @f]
 *
 * The parameters @f$ \rho @f$ and @f$ B @f$ are determined by least-square fitting.
 * Initial values have to be provided.
 * For numerical stability the reflections are split in a high and low resolution set
 * in the fitting procedure. Therefore a resolution cutoff has to be given. Finally
 * the maximum iterations have to be given.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomssf</td><td>&lt;@ref AtomSpecifier atoms to consider for structure_factor&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;file with IAC-to-elementname mapping&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;file with experimental B-factors and occupancies&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution [nm]&gt; </td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup in Hermann-Mauguin format, default: P 1&gt;]</td></tr>
 * <tr><td> \@cif</td><td>&lt;cristallographic information file&gt; </td></tr>
 * <tr><td>[\@density</td><td>&lt;write density to files&gt;]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;convert length unit to Angstrom&gt;]</td></tr>
 * <tr><td>[\@bins</td><td>&lt;number of bins used for fitting&gt;]</td></tr>
 * <tr><td>[\@solvent</td><td>&lt;solvent parameters&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 r_factor
    @topo       ex.top
    @time       0 0.1
    @atomssf    1:CA
    @traj       ex.tr
    @map        ex.map
    @bfactor    ex.bfc
    @resolution 0.1
    @spacegroup P 21 21 21
 #              RVDW RION RSHRINK IACW IRHO IB   RESCUT MAXIT
    @solvent    0.1  0.08 0.11    4    0.34 5.0  0.5    500
    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gio/InIACElementNameMapping.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/gio/InCIF.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/debug.h"

// Additional Clipper Headers
#include "../config.h"
#include "../src/gromos/Exception.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace std;
using namespace gmath;
using namespace bound;

// the data needed
struct fit_data {
  double factor;
  const clipper::CCell * cell;
  double resolution_cutoff;
  const clipper::CHKL_data<clipper::data64::F_phi> * fphi;
  const clipper::CHKL_data<clipper::data64::F_phi> * fphi_solv;
  const vector<CIFData> * cifdata;
};

// this function is for fitting the solvent
// this will work like this:
// k*|Fcalc| - |Fobs| is minimized with k being the normal least-square scaling constant
int Gprime(const gsl_vector * param, void * data, gsl_vector * f) {
  double factor = static_cast<fit_data*>(data)->factor;
  double resolution_cutoff = static_cast<fit_data*>(data)->resolution_cutoff;
  const clipper::CCell & cell = *(static_cast<fit_data*>(data)->cell);
  const clipper::CHKL_data<clipper::data64::F_phi> & fphi = *(static_cast<fit_data*>(data)->fphi);
  const clipper::CHKL_data<clipper::data64::F_phi> & fphi_solv = *(static_cast<fit_data*>(data)->fphi_solv);
  const vector<CIFData> & cifdata = *(static_cast<fit_data*>(data)->cifdata);

  double rho = gsl_vector_get(param, 0);
  double B = gsl_vector_get(param, 1);

  // calculate the scaling constant
  const unsigned int num_bins = 2;
  vector<double> f_calc(cifdata.size());
  vector<double> sum_obs_calc(num_bins, 0.0);
  vector<double> sum_calc_calc(num_bins, 0.0);
  for (unsigned int i = 0; i < cifdata.size(); ++i) {
    clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
    unsigned int bin = 0;
    double inv_res = hkl.invresolsq(cell) * factor * factor;
    double hkl_resolution = sqrt(1.0 / inv_res);
    if (hkl_resolution > resolution_cutoff) bin = 1;

    const double f_obs = cifdata[i].f_obs;

    complex<double> f_calc_model(fphi[hkl]);
    complex<double> f_calc_solv(fphi_solv[hkl]);
    f_calc[i] = abs(f_calc_model + f_calc_solv * rho * exp(-B * inv_res));

    sum_obs_calc[bin] += f_obs * f_calc[i];
    sum_calc_calc[bin] += f_calc[i] * f_calc[i];
  }
  vector<double> k(num_bins, 0.0);
  for (unsigned int bin = 0; bin < num_bins; ++bin) {
    if (sum_calc_calc[bin] != 0.0)
      k[bin] = sum_obs_calc[bin] / sum_calc_calc[bin];
  }

  for (unsigned int i = 0; i < cifdata.size(); ++i) {
    clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
    unsigned int bin = 0;
    double inv_res = hkl.invresolsq(cell) * factor * factor;
    double hkl_resolution = sqrt(1.0 / inv_res);
    if (hkl_resolution > resolution_cutoff) bin = 1;

    const double f_obs = cifdata[i].f_obs;
    gsl_vector_set(f, i, (k[bin] * f_calc[i] - f_obs));
    
  }
  return GSL_SUCCESS;
}
// the jacobian of it i.e.
// k'*|Fcalc| + k*|Fcalc|'
int dGprime_dp(const gsl_vector * param, void * data, gsl_matrix * J) {
  double factor = static_cast<fit_data*>(data)->factor;
  double resolution_cutoff = static_cast<fit_data*>(data)->resolution_cutoff;
  const clipper::CCell & cell = *(static_cast<fit_data*>(data)->cell);
  const clipper::CHKL_data<clipper::data64::F_phi> & fphi = *(static_cast<fit_data*>(data)->fphi);
  const clipper::CHKL_data<clipper::data64::F_phi> & fphi_solv = *(static_cast<fit_data*>(data)->fphi_solv);
  const vector<CIFData> & cifdata = *(static_cast<fit_data*>(data)->cifdata);

  double rho = gsl_vector_get(param, 0);
  double B = gsl_vector_get(param, 1);

  vector<double> f_calc(cifdata.size());
  vector<double> d_fcalc_drho(cifdata.size());
  vector<double> d_fcalc_dB(cifdata.size());
  for (unsigned int i = 0; i < cifdata.size(); ++i) {
    clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
    double inv_res = hkl.invresolsq(cell) * factor * factor;
    complex<double> f_calc_model(fphi[hkl]);
    complex<double> f_calc_solv(fphi_solv[hkl]);
    const double exp_term = exp(-B*inv_res);
    complex<double> cf_calc = f_calc_model + f_calc_model + f_calc_solv*rho*exp_term;
    f_calc[i] = abs(cf_calc);
    const double a_calc = cf_calc.real();
    const double b_calc = cf_calc.imag();
    const double a_solv = f_calc_solv.real();
    const double b_solv = f_calc_solv.imag();

    double term = a_calc * a_solv + b_calc * b_solv;
    d_fcalc_drho[i] = exp_term * term / f_calc[i];
    d_fcalc_dB[i] = -rho * inv_res * d_fcalc_drho[i];
  }

  // calculate the scaling constant
  const unsigned int num_bins = 2;
  vector<double> sum_obs_calc(num_bins, 0.0);
  vector<double> sum_obs_calc_drho(num_bins, 0.0);
  vector<double> sum_obs_calc_dB(num_bins, 0.0);
  vector<double> sum_calc_calc(num_bins, 0.0);
  vector<double> sum_calc_calc_drho(num_bins, 0.0);
  vector<double> sum_calc_calc_dB(num_bins, 0.0);
  for (unsigned int i = 0; i < cifdata.size(); ++i) {
    clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
    unsigned int bin = 0;
    const double inv_res = hkl.invresolsq(cell) * factor * factor;
    const double hkl_resolution = sqrt(1.0 / inv_res);
    if (hkl_resolution > resolution_cutoff) bin = 1;
    const double f_obs = cifdata[i].f_obs;

    sum_obs_calc[bin] += f_obs * f_calc[i];
    sum_obs_calc_drho[bin] += f_obs * d_fcalc_drho[i];
    sum_obs_calc_dB[bin] += f_obs * d_fcalc_dB[i];
    sum_calc_calc[bin] += f_calc[i] * f_calc[i];
    sum_calc_calc_drho[bin] += f_calc[i] * d_fcalc_drho[i];
    sum_calc_calc_dB[bin] += f_calc[i] * d_fcalc_dB[i];
  }
  vector<double> k(num_bins, 0.0);
  vector<double> d_k_drho(num_bins, 0.0);
  vector<double> d_k_dB(num_bins, 0.0);
  for (unsigned int bin = 0; bin < num_bins; ++bin) {
    if (sum_calc_calc[bin] != 0.0) {
      k[bin] = sum_obs_calc[bin] / sum_calc_calc[bin];
      d_k_drho[bin] = (sum_obs_calc_drho[bin] * sum_calc_calc[bin] - 2.0 * sum_obs_calc[bin] * sum_calc_calc_drho[bin]) / (sum_calc_calc[bin] * sum_calc_calc[bin]);
      d_k_dB[bin] = (sum_obs_calc_dB[bin] * sum_calc_calc[bin] - 2.0 * sum_obs_calc[bin] * sum_calc_calc_dB[bin]) / (sum_calc_calc[bin] * sum_calc_calc[bin]);
    }
  }

  for (unsigned int i = 0; i < cifdata.size(); ++i) {
    clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
    unsigned int bin = 0;
    const double inv_res = hkl.invresolsq(cell) * factor * factor;
    const double hkl_resolution = sqrt(1.0 / inv_res);
    if (hkl_resolution > resolution_cutoff) bin = 1;

    gsl_matrix_set(J, i, 0, d_k_drho[bin] * f_calc[i] + k[bin] * d_fcalc_drho[i]);
    gsl_matrix_set(J, i, 1, d_k_dB[bin] * f_calc[i] + k[bin] * d_fcalc_dB[i]);
  }

  return GSL_SUCCESS;
}

int Gprime_and_dGprime_dp(const gsl_vector * param, void * data, gsl_vector * f, gsl_matrix * J) {
  Gprime(param, data, f);
  dGprime_dp(param, data, J);

  return GSL_SUCCESS;
}


int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "map" << "atomssf" << "time" << "bfactor"
          << "resolution" << "spacegroup" << "cif" << "density" << "factor" << "bins" << "solvent";

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
  usage += "\t@cif         <cristallographic information file>\n";
  usage += "\t[@density    <write electron density to special files>]\n";
  usage += "\t[@factor     <convert length unit to Angstrom. default: 10.0>]\n";
  usage += "\t[@bins       <number of bins used for Fobs/Fcalc fitting>]\n";
  usage += "\t[@solvent    <solvent parameters. Many: see help>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);

  try {
    Arguments args(argc, argv, knowns, usage);
    double factor = args.getValue<double>("factor", false, 10.0);
    unsigned int num_bins = args.getValue<int>("bins", false, 1);

    // Hardcoded B-factor conversion factor.
    const double sqpi2=(M_PI*M_PI*8.0);

    // Get Spacegroup Data or default to no symmetry (P 1)
    string spgrdata("P 1");
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
    auto_ptr<clipper::Spgr_descr> spgrinit;
    try {
     spgrinit = auto_ptr<clipper::Spgr_descr>(new clipper::Spgr_descr(spgrdata, clipper::Spgr_descr::HM));
    } catch(clipper::Message_fatal & msg) {
      throw gromos::Exception(argv[0], "Invalid spacegroup: " + msg.text());
    }
    clipper::CSpacegroup spgr(clipper::String("base spgr"), clipper::Spacegroup(*spgrinit));

    // Get resolution as a double
    double resolution, resolution_max = std::numeric_limits<double>::max();
    {
      Arguments::const_iterator iter = args.lower_bound("resolution");
      Arguments::const_iterator to = args.upper_bound("resolution");
      if (iter == to) {
        throw gromos::Exception(argv[0], "Please give at least one resolution parameter");
      }
      istringstream iss (iter->second);
      if (!(iss >> resolution)) {
        throw gromos::Exception(argv[0],
                "Resolution parameter not numeric");
      }
      ++iter;
      istringstream iss2 (iter->second);
      if (iter != to) {
        if (!(iss2 >> resolution_max)) {
          throw gromos::Exception(argv[0],
                  "Resolution parameter not numeric");
        }
      }
    }

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    GromosForceField gff(it.forceField());
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
    map<int, string> gacmapping = mapfile.getData();

    // Get experimental Bfactors and occupancy
    vector<BFactorOccupancyData> bfoc;
    bool has_bfactor = false;
    if (args.count("bfactor") == 1) {
      InBFactorOccupancy bfac_file(args["bfactor"]);
      bfoc = bfac_file.getData();
      has_bfactor = true;
    }

    // Read in cif-file for generating the structure factor list (using clipper)
    InCIF ciffile(args["cif"]);
    vector<CIFData> init_cifdata(ciffile.getData());

    bool write_density = false;
    if (args.count("density") >= 0) {
      write_density = true;
    }

    bool solvent_corr = false;
    double vdw_probe = 0.14;
    double ion_probe = 0.08;
    double r_shrink = 0.11;
    int water_iac = 4;
    double rho_init = 0.34;
    double B_init = 5.0;
    double resolution_cutoff = 0.5;
    unsigned int max_iter = 500;
    {
      Arguments::const_iterator it = args.lower_bound("solvent"),
              to = args.upper_bound("solvent");
      if (it != to) {
        solvent_corr = true;
        if (!(istringstream(it->second) >> vdw_probe))
          throw gromos::Exception(argv[0], "vdW probe radius not numeric.");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> ion_probe))
          throw gromos::Exception(argv[0], "ion probe radius not numeric.");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> r_shrink))
          throw gromos::Exception(argv[0], "shrink radius not numeric");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> water_iac))
          throw gromos::Exception(argv[0], "water IAC not numeric.");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> rho_init))
          throw gromos::Exception(argv[0], "initial rho value not numeric");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> B_init))
          throw gromos::Exception(argv[0], "initial B value not numeric");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> resolution_cutoff))
          throw gromos::Exception(argv[0], "resolution cutoff not numeric");
        if (++it == to)
          throw gromos::Exception(argv[0], "not enough solvent arguments");
        if (!(istringstream(it->second) >> max_iter))
          throw gromos::Exception(argv[0], "maximum iterations not numeric");
      }
    }

    //===========================
    // loop over all trajectories
    InG96 ic;

    cout << "#" << setw(14) << "time" << setw(8) << "R_tot" << setw(8) << "R[0]" << setw(8) << "k[0]" << " ..." << endl;
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

        // create the cell
        const double a = sys.box().K().abs();
        const double b = sys.box().L().abs();
        const double c = sys.box().M().abs();
        DEBUG(5, "Cell: " << a << " " << b << " " << c)
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
        clipper::CHKL_data<clipper::data64::F_phi> fphi_solv(hkls);
        clipper::CHKL_data<clipper::data64::F_phi> fphi_print(hkls);
        vector<CIFData> cifdata;
        for (unsigned int i = 0; i < init_cifdata.size(); ++i) {
          clipper::HKL hkl(init_cifdata[i].h, init_cifdata[i].k, init_cifdata[i].l);
          const double inv_res = hkl.invresolsq(cell) * factor * factor;
          const double hkl_resolution = sqrt(1.0 / inv_res);
          if (hkl_resolution >  resolution && hkl_resolution < resolution_max) {
            cifdata.push_back(init_cifdata[i]);
          }
        }


        // Fill Clipper Atom list
        // we do this insight the loop due to solvent molecules!
        vector<clipper::Atom> atomvec;
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

        if (solvent_corr) {
          const clipper::Grid_sampling solvent_grid(spgr, cell, reso, 1.5);
          clipper::Xmap<clipper::ftype64> mask(spgr, cell, solvent_grid);
          mask = 1.0;
          set<int> mask_indices;
          for (int i = 0; i < calcatoms.size(); i++) {
            double vdw_radius = 0.0;
            LJType lj(gff.ljType(AtomPair(calcatoms.iac(i), water_iac)));
            if (lj.c6() >= 1.0E-20) {
              vdw_radius = exp(log(2.0 * lj.c12() / lj.c6()) / 6.0);
            }

            vdw_radius /= 2.0;

            if (vdw_radius && calcatoms.charge(i) != 0.0)
              vdw_radius += ion_probe;
            else
              vdw_radius += vdw_probe;

            double cutoff2 = vdw_radius * vdw_radius * factor * factor;

            clipper::Grid_range gd(cell, solvent_grid, vdw_radius * factor);
            clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(cell);
            clipper::Coord_grid g0 = uvw.coord_grid(solvent_grid) + gd.min();
            clipper::Coord_grid g1 = uvw.coord_grid(solvent_grid) + gd.max();
            // loop over atom's grid
            clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
            i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(mask, g0);
            for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
              for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
                for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                  const clipper::Coord_orth dr(atoms[i].coord_orth()-iw.coord_orth());
                  if (dr.lengthsq() < cutoff2) {
                    mask_indices.insert(iw.index());
                  }
                }
              }
            } // loop over grid
          } // loop over atoms
          // shrink the mask
          double cutoff_shrink = r_shrink * r_shrink * factor * factor;
          for (set<int>::const_iterator it = mask_indices.begin(), to = mask_indices.end();
                  it != to; ++it) {
            mask.set_data(*it, 0.0);
            clipper::Coord_grid grid_point(mask.coord_of(*it));
            clipper::Coord_orth grid_point_orth(grid_point.coord_frac(solvent_grid).coord_orth(cell));
            clipper::Grid_range gd(cell, solvent_grid, r_shrink * factor);
            clipper::Coord_grid g0 = mask.coord_of(*it) + gd.min();
            clipper::Coord_grid g1 = mask.coord_of(*it) + gd.max();
            // loop over atom's grid
            clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
            i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(mask, g0);
            for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
              for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
                for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
                  // check for cell cutoff
                  if ((iw.coord_orth() - grid_point_orth).lengthsq() < cutoff_shrink &&
                      mask_indices.find(iw.index()) == mask_indices.end()) {
                    mask.set_data(*it, 1.0);
                    goto out;
                  }
                }
              }
            } // loop over grid
            out: ;
          } // loop over mask

          mask.fft_to(fphi_solv);
          fit_data fd;
          fd.cell = &cell;
          fd.cifdata = &cifdata;
          fd.fphi = &fphi;
          fd.fphi_solv = &fphi_solv;
          fd.resolution_cutoff = resolution_cutoff;
          fd.factor = factor;

          gsl_multifit_function_fdf f;
          f.f = &Gprime;
          f.df = &dGprime_dp;
          f.fdf = &Gprime_and_dGprime_dp;
          f.n = cifdata.size();
          f.p = 2;
          f.params = &fd;

          double param_init[2] = {rho_init, B_init};
          gsl_vector_view param = gsl_vector_view_array(param_init, 2);

          const gsl_multifit_fdfsolver_type * solver_type = gsl_multifit_fdfsolver_lmsder;
          gsl_multifit_fdfsolver * solver = gsl_multifit_fdfsolver_alloc(solver_type, cifdata.size(), 2);
          gsl_multifit_fdfsolver_set(solver, &f, &param.vector);

          unsigned int iter = 0;
          int status;
          double rho, B;
          do {
            iter++;
            status = gsl_multifit_fdfsolver_iterate(solver);
            if (status) break;
            status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-4, 1e-4);
            rho = gsl_vector_get(solver->x, 0);
            B = gsl_vector_get(solver->x, 1);
          } while(status == GSL_CONTINUE && iter < max_iter);
          rho = gsl_vector_get(solver->x, 0);
          B = gsl_vector_get(solver->x, 1);
          cout << "# rho is " << rho << ", B is " << B << " after " << iter << " iterations." << endl;
          gsl_multifit_fdfsolver_free(solver);

          for (unsigned int i = 0; i < cifdata.size(); ++i) {
            clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
            double inv_res = hkl.invresolsq(cell) * factor * factor;
            complex<double> f_calc_model(fphi[hkl]);
            complex<double> f_calc_solv(fphi_solv[hkl]);
            complex<double> f_calc = f_calc_model + f_calc_solv * rho * exp(-B * inv_res);
            fphi.set_data(hkl, f_calc);
          }
        }

        if (write_density)
          fphi_print = complex<double>(0.0,0.0);

        double reso_min = numeric_limits<double>::max(), reso_max = numeric_limits<double>::min();
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
          double hkl_resolution = sqrt(1.0 / hkl.invresolsq(cell)) / factor;
          reso_min = min(reso_min, hkl_resolution);
          reso_max = max(reso_max, hkl_resolution);
        }

        // calculate the scaling constant
        vector<double> sum_obs_calc(num_bins, 0.0);
        vector<double> sum_calc_calc(num_bins, 0.0);
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
          unsigned int bin = 0;
          if (num_bins > 1) {
            double hkl_resolution = sqrt(1.0 / hkl.invresolsq(cell)) / factor;
            bin = (hkl_resolution - reso_min) * num_bins / (reso_max - reso_min);
          }
          const double f_obs = cifdata[i].f_obs;
          const double f_calc = fphi[hkl].f();
          sum_obs_calc[bin] += f_obs * f_calc;
          sum_calc_calc[bin] += f_calc * f_calc;
        }
        vector<double> k(num_bins, 0.0);
        for(unsigned int bin = 0; bin < num_bins; ++bin) {
          if (sum_calc_calc[bin] != 0.0)
            k[bin] = sum_obs_calc[bin] / sum_calc_calc[bin];
        }
        // and calculate R
        double sum_dev_f = 0.0;
        double sum_obs = 0.0;
        vector<double> sum_dev_f_bin(num_bins, 0.0);
        vector<double> sum_obs_bin(num_bins, 0.0);
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
          unsigned int bin = 0;
          if (num_bins > 1) {
            double hkl_resolution = sqrt(1.0 / hkl.invresolsq(cell)) / factor;
            bin = (hkl_resolution - reso_min) * num_bins / (reso_max - reso_min);
          }
          const double f_obs = cifdata[i].f_obs;
          const double f_calc = fphi[hkl].f();
          const double dev = fabs(f_obs - k[bin]*f_calc);

          sum_dev_f += dev;
          sum_dev_f_bin[bin] += dev;
          sum_obs += f_obs;
          sum_obs_bin[bin] += f_obs;
          if (write_density)
            fphi_print.set_data(hkl, clipper::data64::F_phi(2.0 * f_obs - k[bin]*f_calc, fphi[hkl].phi()));
        }

        const double R = sum_dev_f / sum_obs;
        vector<double> R_bin(num_bins, 0.0);
        for(unsigned int bin = 0; bin < num_bins; ++bin) {
          if (sum_obs_bin[bin] != 0.0)
            R_bin[bin] = sum_dev_f_bin[bin] / sum_obs_bin[bin];
        }

        cout.precision(8);
        cout << time;
        cout.precision(5);
        cout << setw(8) << R;
        for(unsigned int bin = 0; bin < num_bins; ++bin) {
          cout << setw(8) << R_bin[bin] << setw(8) << k[bin];
        }
        cout << endl;

        // write the electron density
        if (write_density) {
          {
            ostringstream file_name;
            file_name << "FRAME_DENSITY_" << setfill('0') << setw(6) << frame_number << ".ccp4";
            const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
            clipper::Xmap<clipper::ftype64> density(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);
            density.fft_from(fphi_print);

            for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = density.first();
                    !ix.last(); ix.next())
              density[ix] /= density.multiplicity(ix.coord());

            clipper::CCP4MAPfile mapfile;

            mapfile.open_write(file_name.str());
            mapfile.export_xmap(density);
            mapfile.close_write();
          }
          {
            ostringstream file_name;
            file_name << "FRAME_CALCDENSITY_" << setfill('0') << setw(6) << frame_number << ".ccp4";
            const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
            clipper::Xmap<clipper::ftype64> density(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);
            density.fft_from(fphi);

            for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = density.first();
                    !ix.last(); ix.next())
              density[ix] /= density.multiplicity(ix.coord());

            clipper::CCP4MAPfile mapfile;

            mapfile.open_write(file_name.str());
            mapfile.export_xmap(density);
            mapfile.close_write();
          }
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
