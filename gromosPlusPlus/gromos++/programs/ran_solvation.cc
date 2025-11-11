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
 * @file ran_solvation.cc
 * Solvates a solute by randomly placing solvent molecules around it
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ran_solvation
 * @section ran_solvation Solvates a solute by randomly placing solvent molecules around it
 * @author @ref dt
 * @date 29-11-2004
 *
 * When simulating a molecule in solution or in a crystal containing solvent 
 * molecules, the atomic coordinates of the solvent molecules are to be 
 * generated, if they are not available from experiment. Alternatively to 
 * @ref sim_box (which solvates a solute in a pre-equilibrated box of a 
 * molecular liquid, see section V-2.12) ran_solvation can solvate a solute in 
 * a mixture consisting of an unlimited number of components. The program 
 * places the solute in the center of a box (rectangular or truncated 
 * octahedron) and generates a random distribution of solvent molecules around
 * it, which are placed in a random orientation. The total number of solvent
 * molecules is calculated based on the specified solvent density and the molar
 * fractions. When calculating the solvent density, the excluded volume of the
 * solute is taken into account as specified by the user. Ran_solvation checks
 * separately for solute-solvent and solvent-solvent overlap after every
 * insertion. If any solute-solvent or solvent-solvent interatomic distance is
 * smaller than the respective threshold distance as specified by the user, the
 * insertion trial is rejected. Note that for the optional input flags, either 
 * the box-size or the minimum solute-to-wall distance should be specified.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo_u</td><td>&lt;molecular topology file for the solute&gt; </td></tr>
 * <tr><td> \@pos_u</td><td>&lt;input coordinate file for the solute&gt; </td></tr>
 * <tr><td> \@sev</td><td>&lt;solvent-excluded volume of solute (nm^3)&gt; </td></tr>
 * <tr><td> \@topo_v</td><td>&lt;list of (single molecule) molecular topology files of solvents: topo_solv1 topo_solv2 ...&gt; </td></tr>
 * <tr><td> \@pos_v</td><td>&lt;list of (single molecule) coordinate files of solvents: pos_solv1 pos_solv2 ...&gt; </td></tr>
 * <tr><td> \@molf_v</td><td>&lt;mole fraction of each solvent: molf_solv1 molf_solv2 ...&gt; </td></tr>
 * <tr><td> \@dens</td><td>&lt;mass density of solvent mixture (kg/m^3)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary condiation of the resulting box&gt; </td></tr>
 * <tr><td> [\@minwall</td><td>&lt;minimum solue-to-wall distance(s)&gt;] </td></tr>
 * <tr><td> [\@boxsize</td><td>&lt;length of box-edge(s)&gt;] </td></tr>
 * <tr><td> [\@thresh_u</td><td>&lt;threshold distance in overlap check (solute - solvent) ; default: 0.40 nm&gt;] </td></tr>
 * <tr><td> [\@thresh_v</td><td>&lt;threshold distance in overlap check (solvent - solvent) ; default: 0.20 nm&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ran_solvation
    @topo_u    lyso_ph2mta53a6.top
    @pos_u     lyso_ph2sx_4.g96
    @sev       16
    @topo_v    ic4.top urea.top h2o.top
    @pos_v     ic4.g96 urea.g96 h2o.g96
    @molf_v    0.10    0.15     0.75
    @dens      900
    @pbc       t
    @minwall   1.5
    @thresh_u  0.4
    @thresh_v  0.22
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <unistd.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InTopology.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace gmath;
using namespace bound;
using namespace utils;


// defining some functions
void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as);
double calc_max_size(System &sys, AtomSpecifier &as, int dim);
double calc_mol_radius(System &sys);
double calc_mass_solv(const vector<string> &tops, const int i);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo_u" << "pos_u" << "sev" << "topo_v" << "pos_v" << "dens"
          << "molf_v" << "pbc" << "minwall" << "boxsize" << "thresh_u"
          << "thresh_v";

  string usage = "# " + string(argv[0]);
  usage += "\n";
  usage += "\t@topo_u      <molecular topology file for the solute>\n";
  usage += "\t@pos_u       <input coordinate file for the solute>\n";
  usage += "\t@sev         <solvent-excluded volume of solute (nm^3)>\n";
  usage += "\n";
  usage += "\t@topo_v      <list of (single molecule) molecular topology files of solvents: topo_solv1 topo_solv2 ...>\n";
  usage += "\t@pos_v       <list of (single molecule) coordinate files of solvents: pos_solv1  pos_solv2 ...>\n";
  usage += "\t@molf_v      <mole fraction of each solvent: molf_solv1 molf_solv2 ...>\n";
  usage += "\t@dens        <mass density of solvent mixture (kg/m^3)>\n";
  usage += "\n";
  usage += "\t@pbc         <periodic boundary condiation of the resulting box>\n";
  usage += "\t@minwall     <minimum solue-to-wall distance(s)>\n";
  usage += "\t@boxsize     <length of box-edge(s)>\n";
  usage += "\n";
  usage += "\t@thresh_u    <threshold distance in overlap check (solute  - solvent) ; default: 0.40 nm>\n";
  usage += "\t@thresh_v    <threshold distance in overlap check (solvent - solvent) ; default: 0.20 nm>\n";

  srand(time(NULL));

  try {


    Arguments args(argc, argv, knowns, usage);
    Arguments::const_iterator iter;

    // topo_pro:
    string topo_pro = args.getValue<std::string > ("topo_u", true, "");
    // insx_pro:
    string insx_pro = args.getValue<std::string > ("pos_u", true, "");
    // sev:
    double sev = args.getValue<double>("sev");

    // consistency check of solvent specification:
    if (args.count("topo_v") != args.count("pos_v")) {
      std::ostringstream os;
      os << "\nFor every solvent, you need to specify a topology and a "
              << "coordinate file\n"
              << "You have specified " << args.count("topo_v") << " topolog";
      if (args.count("topo_v") == 1) os << "y";
      else os << "ies";
      os << " and "
              << args.count("pos_v") << " coordinate file";
      if (args.count("pos_v") != 1) os << "s";
      os << "\n";
      throw gromos::Exception("ran_solvation", os.str());
    }

    // topo_solv:
    args.check("topo_v", 1);
    vector<string> tops;
    iter = args.lower_bound("topo_v");
    while (iter != args.upper_bound("topo_v")) {
      tops.push_back(iter->second.c_str());
      ++iter;
    }
    // insx_solv:
    args.check("pos_v", 1);
    vector<string> insxs;
    iter = args.lower_bound("pos_v");
    while (iter != args.upper_bound("pos_v")) {
      insxs.push_back(iter->second.c_str());
      ++iter;
    }
    // mol_frac:
    args.check("molf_v", 1);
    vector<double> molfsv;
    iter = args.lower_bound("molf_v");
    while (iter != args.upper_bound("molf_v")) {
      molfsv.push_back(atof(iter->second.c_str()));
      ++iter;
    }
    // densit:
    double densit = args.getValue<double>("dens", true);

    // minwall:    
    // read the minimum solute to wall distances.
    // three possibilities:
    // 1. nothing specified: the user will specify the box size
    //    using the @boxsize flag
    // 2. one number specified: the box will be cubic (pbc: r) or a
    //    trunc. oct. (pbc: t) with the size defined by maximum solute
    //    atom-atom distance plus twice this number
    // 3. three numbers specified: the box will be rectangular. The
    //    solute is rotated and the specified numbers are added to the
    //    maximum distances in the resulting dimensions
    vector<double> minwall;
    //int calculate_dimensions = 0;
    {
      Arguments::const_iterator iter = args.lower_bound("minwall"),
              to = args.upper_bound("minwall");
      while (iter != to) {
        minwall.push_back(atof(iter->second.c_str()));
        ++iter;
      }
      //if (minwall.size()) calculate_dimensions = 1;
    }

    // read the boxsize
    // three posibilities:
    // 1. nothing specified: the box size is determined from solute 
    //    geometry and minwall specifications
    // 2. one number specified: the box will be cubic (pbc: r) or a trunc.
    //    oct. (pbc:t)
    // 3. three numbers specified: the box will be rectangular
    vector<double> boxsize;
    {
      Arguments::const_iterator iter = args.lower_bound("boxsize"),
              to = args.upper_bound("boxsize");
      while (iter != to) {
        boxsize.push_back(atof(iter->second.c_str()));
        ++iter;
      }
      if (boxsize.size() == 1 || boxsize.size() == 2)
        boxsize.resize(3, boxsize[0]);
    }
    if (boxsize.size() && minwall.size())
      throw (gromos::Exception("ran_solvation",
            "Cannot specify both boxsize and minwall."));
    // pbc:
    InTopology it(args["topo_u"]);
    System solu(it.system());
    Box box;
    Boundary *pbc = BoundaryParser::boundary(solu, args);
    //int truncoct = 0, rectbox = 0, cubic = 0;
    double size_corr = 1.0;
    double fac_vol = 1.0;
    switch (pbc->type()) {
      case('t'):
        //truncoct = 1;
        //rectbox = 0;
        size_corr = 2.0 * sqrt(3.0) / 3.0;
        fac_vol = 0.5;
        if (minwall.size() > 1)
          throw (gromos::Exception("ran_solvation",
                "For truncated octahedral boxes you can only specify one number for @minwall"));
        if (minwall.size() == 0)
          if (boxsize[0] != boxsize[1] || boxsize[0] != boxsize[2])
            throw (gromos::Exception("ran_solvation",
                 "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));
        box.setNtb(gcore::Box::truncoct);
        break;
      case('r'):
        //truncoct = 0;
        //rectbox = 1;
        fac_vol = 1.0;
        //if (minwall.size() == 1) cubic = 1;
        //else if (minwall.size() == 0 &&
        //        (boxsize[0] == boxsize[1] && boxsize[0] == boxsize[2]))
          //cubic = 1;
        box.setNtb(gcore::Box::rectangular);
        break;
      case('v'):
        throw (gromos::Exception("ran_solvation",
                "Why are you running this program if @pbc is vacuum?"));
        break;
      default:
        stringstream msg;
        msg << "Periodic boundary condition " << pbc->type() << " is not supported.";
        throw gromos::Exception("ran_solvation", msg.str());
        break;
    }

    // reading the threshold values
    //
    // thresh_pro    
    iter = args.lower_bound("thresh_u");
    double thresh_pro = (iter != args.upper_bound("thresh_u")) ? atof(iter->second.c_str())
            : 0.4;
    thresh_pro *= thresh_pro;
    // thresh_solv
    iter = args.lower_bound("thresh_v");
    double thresh_solv = (iter != args.upper_bound("thresh_v")) ? atof(iter->second.c_str())
            : 0.20;
    thresh_solv *= thresh_solv;

    // Compute vol_box:
    //  
    InG96 ic;
    ic.open(args["pos_u"]);
    ic >> solu;
    ic.close();

    solu.box() = box;

    // shift molecule to cog
    double radius_pro = calc_mol_radius(solu);
    Vec shiftcog = fit::PositionUtils::shiftToCog(&solu);

    // there is only need to rotate if minwall has three elements
    // store the maximum dimensions in max_dim
    vector<double> max_dim;
    AtomSpecifier as(solu);
    if (minwall.size() == 3) {
      rotate_solute(solu, max_dim, as);
      // calculate the box dimensions
      for (int i = 0; i < 3; i++)
        boxsize.push_back(max_dim[i] + 2 * minwall[i]);
    } else if (boxsize.size() == 0) {
      max_dim.push_back(calc_max_size(solu, as, 3));
      // calculate box dimensions
      for (int i = 0; i < 3; i++)
        boxsize.push_back(size_corr * (max_dim[0] + 2 * minwall[0]));
    }

    solu.box().K()[0] = boxsize[0];
    solu.box().L()[1] = boxsize[1];
    solu.box().M()[2] = boxsize[2];

    Vec box_mid(solu.box().K()[0] / 2.0, solu.box().L()[1] / 2.0, solu.box().M()[2] / 2.0);
    fit::PositionUtils::translate(&solu, box_mid);
    double vol_cell = fac_vol * boxsize[0] * boxsize[1] * boxsize[2];

    // Compute vol_solv:
    double vol_solv = vol_cell - sev;

    // computing the masses of the solvents and
    vector<double> mass_solv;
    for (unsigned int i = 0; i < tops.size(); i++)
      mass_solv.push_back(calc_mass_solv(tops, i));

    // Computing the number of molecules of each solvent required to match densit:
    double tmp = 0;
    for (unsigned int i = 0; i < tops.size(); i++)
      tmp += molfsv[i] * mass_solv[i];
    vector<double> nr_solv; // number of molecules for each type of slv
    vector<double> tot_mass_solv;
    for (unsigned int i = 0; i < tops.size(); i++) {
      nr_solv.push_back(int(densit * vol_solv * molfsv[i] / tmp));
      tot_mass_solv.push_back(mass_solv[i] * nr_solv[i]);
    }

    //
    // report computational specs:
    // 
    cerr << "Cell volume: " << vol_cell << endl
            << "SEV of solute: " << sev << endl
            << "Solvent volume: " << vol_solv << endl
            << "Input solvent density: " << densit << endl;
    for (unsigned int i = 0; i < tops.size(); i++)
      cerr << "# solvent " << i + 1 << ": " << nr_solv[i] << endl
            << "Total mass of solvent " << i + 1 << ": " << tot_mass_solv[i] << endl;
    cerr << endl << "CELL DIMENSIONS: " << endl
            << "X: " << boxsize[0]
            << "\t\tY: " << boxsize[1]
            << "\t\tZ: " << boxsize[2] << endl
            << "PBC: " << args["pbc"] << endl << endl;

    cerr << "Now, sleeping for 10 seconds... "
            << "there is still time for a ctrl-C!" << endl;
    sleep(10);

    //
    // MAIN PROCESS:
    //
    // loop over the number of solvent topologies.
    for (unsigned int tcnt = 0; tcnt < tops.size(); tcnt++) {

      //read topologies again
      InTopology it(tops[tcnt]);
      System smol(it.system());

      // read single molecule coordinates...
      InG96 ic;
      ic.open(insxs[tcnt]);
      ic >> smol;
      ic.close();

      // single molecule coordinates are translated to reference frame 
      // with origin at centre of mass
      fit::PositionUtils::shiftToCog(&smol);

      //loop over the number of desired mol
      double nr_mol = nr_solv[tcnt];
      for (unsigned int i = 0; i < nr_mol; i++) {
        solu.addMolecule(smol.mol(0));

UGLY_GOTO:
        // trial position of solvent molecule
        //
        Vec rpos;
        for (int d = 0; d < 3; d++) {
          int r = rand();
          rpos[d] = double(solu.box().K()[0] * r) / double(RAND_MAX);
        }
        // if outside box then discard
        if (!pbc->inBox(box_mid, rpos, solu.box())) {
          goto UGLY_GOTO;
        }

        // trial orientation of solvent molecule
        //
        int r = rand();
        int r_axis = int(double(3.0 * r) / (double(RAND_MAX) + 1)) + 1; // i.e. [1,2,3]
        double phi = 2.0 * gmath::physConst.get_pi() * double(r) / (double(RAND_MAX) + 1);
        double cosp = cos(phi), sinp = sin(phi);

        // rotate and put the molecule in the box
        if (r_axis == 1)
          for (int j = 0; j < smol.mol(0).numAtoms(); j++) {
            solu.mol(solu.numMolecules() - 1).pos(j)[0] = smol.mol(0).pos(j)[0] + rpos[0];
            solu.mol(solu.numMolecules() - 1).pos(j)[1] = (smol.mol(0).pos(j)[1] * cosp -
                    smol.mol(0).pos(j)[2] * sinp) + rpos[1];
            solu.mol(solu.numMolecules() - 1).pos(j)[2] = (smol.mol(0).pos(j)[2] * cosp +
                    smol.mol(0).pos(j)[1] * sinp) + rpos[2];
          }
        else if (r_axis == 2)
          for (int j = 0; j < smol.mol(0).numAtoms(); j++) {
            solu.mol(solu.numMolecules() - 1).pos(j)[0] = (smol.mol(0).pos(j)[0] * cosp +
                    smol.mol(0).pos(j)[2] * sinp) + rpos[0];
            solu.mol(solu.numMolecules() - 1).pos(j)[1] = smol.mol(0).pos(j)[1] + rpos[1];
            solu.mol(solu.numMolecules() - 1).pos(j)[2] = (smol.mol(0).pos(j)[2] * cosp -
                    smol.mol(0).pos(j)[0] * sinp) + rpos[2];
          } else if (r_axis == 3)
          for (int j = 0; j < smol.mol(0).numAtoms(); j++) {
            solu.mol(solu.numMolecules() - 1).pos(j)[0] = (smol.mol(0).pos(j)[0] * cosp -
                    smol.mol(0).pos(j)[1] * sinp) + rpos[0];
            solu.mol(solu.numMolecules() - 1).pos(j)[1] = (smol.mol(0).pos(j)[1] * cosp +
                    smol.mol(0).pos(j)[0] * sinp) + rpos[1];
            solu.mol(solu.numMolecules() - 1).pos(j)[2] = smol.mol(0).pos(j)[2] + rpos[2];
          }

        //checking overlap:
        // first with main solute (protein)
        radius_pro += thresh_pro;
        radius_pro *= radius_pro;
        if ((box_mid - rpos).abs2() <= radius_pro) {
          for (int l = 0; l < solu.mol(0).numAtoms(); l++) {
            for (int m = 0; m < smol.mol(0).numAtoms(); m++) {
              // no pbc check because protein certainly within box
              if ((solu.mol(0).pos(l) - solu.mol(solu.numMolecules() - 1).pos(m)).abs2() < thresh_pro) {
                goto UGLY_GOTO;
              }

            }
          }
        }//if rpos in sphere around protein

        // then with all other solvent molecules
        for (int k = 1; k < solu.numMolecules() - 1; k++) {
          for (int l = 0; l < solu.mol(k).numAtoms(); l++) {
            for (int m = 0; m < smol.mol(0).numAtoms(); m++) {
              if ((solu.mol(k).pos(l) -
                      (pbc->nearestImage(solu.mol(k).pos(l), solu.mol(solu.numMolecules() - 1).pos(m), solu.box()))).abs2()
                      < thresh_solv) {
                goto UGLY_GOTO;
              }
            }
          }
        }

        cerr << (i + 1) << " of " << nr_mol << " copies of solvent molecule of type " << tcnt + 1
                << " already in the box. (Total number of molecules = " << solu.numMolecules() << ")." << endl;
      }
      cerr << "Box now with: " << solu.numMolecules() << " molecules" << endl;
    }

    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << solu;

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


// function bodies

void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as) {
  max.resize(3);

  // this will be a two stage function.
  // First calculate the maximum distance between any two solute atoms
  max[2] = calc_max_size(sys, as, 3);

  // rotate the solute such that the atoms in as are along the z-axis
  Vec v = sys.mol(as.mol(0)).pos(as.atom(0))
          - sys.mol(as.mol(1)).pos(as.atom(1));
  double r = v.abs();
  double r_yz = sqrt(v[1] * v[1] + v[2] * v[2]);

  // the rototation matrix is the product of two rotations
  // 1. around the x-axis by theta
  //    with sin(theta) = v[2]/r_yz; cos(theta) = align[2]/r_yz
  // 2. around the y-axis by phi
  // with sin(phi) = align[0]/r; cos(phi) = r_yz / r
  Matrix rot1(Vec(r_yz / r, 0, v[0] / r),
          Vec(-v[0] * v[1] / r / r_yz, v[2] / r_yz, v[1] / r),
          Vec(-v[0] * v[2] / r / r_yz, -v[1] / r_yz, v[2] / r));

  for (int m = 0; m < sys.numMolecules(); m++)
    for (int a = 0; a < sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot1 * sys.mol(m).pos(a);
  // and take along any solvent
  for (int i = 0; i < sys.sol(0).numPos(); i++) {
    sys.sol(0).pos(i) = rot1 * sys.sol(0).pos(i);
  }

  // calculate the maximum distance in the x-y-plane
  max[1] = calc_max_size(sys, as, 2);

  // rotate the solute around the z-axis, such that the atoms in as are
  // along the y-axis, this is done by a rotation around psi with
  // sin(psi) = x/r_xy; cos(psi) = y/r_xy;
  v = sys.mol(as.mol(2)).pos(as.atom(2)) - sys.mol(as.mol(3)).pos(as.atom(3));
  double r_xy = sqrt(v[0] * v[0] + v[1] * v[1]);

  Matrix rot2(Vec(+v[1] / r_xy, v[0] / r_xy, 0),
          Vec(-v[0] / r_xy, v[1] / r_xy, 0),
          Vec(0, 0, 1));

  for (int m = 0; m < sys.numMolecules(); m++)
    for (int a = 0; a < sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot2 * sys.mol(m).pos(a);

  // and take along any solvent
  for (int i = 0; i < sys.sol(0).numPos(); i++) {
    sys.sol(0).pos(i) = rot2 * sys.sol(0).pos(i);
  }

  // finally we calculate the maximum distance in the x-direction
  max[0] = calc_max_size(sys, as, 1);
}

double calc_max_size(System &sys, AtomSpecifier &as, int dim) {
  // calculate the longest distance between solute atoms considering
  // the first dim dimensions.
  double max2 = 0.0;
  double d2 = 0;
  int max_m1 = 0, max_m2 = 0, max_a1 = 0, max_a2 = 0;

  for (int m1 = 0; m1 < sys.numMolecules(); m1++) {
    for (int a1 = 0; a1 < sys.mol(m1).numAtoms(); a1++) {
      Vec current = sys.mol(m1).pos(a1);
      for (int m2 = m1; m2 < sys.numMolecules(); m2++) {
        int start = 0;
        if (m1 == m2) start = a1;
        for (int a2 = start; a2 < sys.mol(m2).numAtoms(); a2++) {
          d2 = 0.0;
          for (int i = 0; i < dim; i++) {
            d2 += (current[i] - sys.mol(m2).pos(a2)[i]) *
                    (current[i] - sys.mol(m2).pos(a2)[i]);
          }
          if (d2 > max2) {
            max2 = d2;
            max_m1 = m1;
            max_m2 = m2;
            max_a1 = a1;
            max_a2 = a2;
          }
        }
      }
    }
  }
  as.addAtomStrict(max_m1, max_a1);
  as.addAtomStrict(max_m2, max_a2);

  return sqrt(max2);
}

double calc_mol_radius(System &sys) {
  // calculate the radius of the solute molecule
  double d2 = 0;

  Vec cog(0.0, 0.0, 0.0);
  int counter = 0;
  for (int m2 = 0; m2 < sys.numMolecules(); m2++) {

    for (int a2 = 0; a2 < sys.mol(m2).numAtoms(); a2++) {
      cog += sys.mol(m2).pos(a2);
    }
    counter += sys.mol(m2).numAtoms();
  }
  cog /= counter;
  for (int m2 = 0; m2 < sys.numMolecules(); m2++) {

    for (int a2 = 0; a2 < sys.mol(m2).numAtoms(); a2++) {
      double t = (cog - sys.mol(m2).pos(a2)).abs2();

      if (d2 < t) d2 = t;
    }
  }
  return sqrt(d2);
}

double calc_mass_solv(const vector<string> &tops, const int i) {
    
  // setting some constants
  const double nano_i = 1 / gmath::physConst.get_nano();
  const double fac_amu2kg = gmath::physConst.get_atomic_mass_unit() * nano_i * nano_i * nano_i ;

  // calculates the total mass of type i solvent
  InTopology itA(tops[i]);
  System solv_A(itA.system());

  double mass_A = 0;
  if (solv_A.numMolecules() != 1) {
    throw (gromos::Exception("ran_solvation",
            "Topology of first solvent has more than one molecule in SOLUTEATOMS block."));
  } else {
    for (int j = 0; j < solv_A.mol(0).numAtoms(); j++) {
      mass_A += solv_A.mol(0).topology().atom(j).mass();
    }
    mass_A *= fac_amu2kg;
  }
  return mass_A;
}


