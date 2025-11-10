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
 * @file ran_box.cc
 * Create a condensed phase system of any composition
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ran_box
 * @section ran_box Create a condensed phase system of any composition
 * @author @ref dt
 * @date 7-6-07
 *
 * When simulating a molecular liquid, a starting configuration for the solvent
 * molecules has to be generated. Program ran_box generates a starting
 * configuration for the simulation of mixtures consisting of an unlimited 
 * number of components. The molecules are randomly placed in a cubic or a
 * truncated octahedron box, in a random orientation. Note that for the
 * generation of a starting configuration for the simulation of pure liquids
 * and binary mixtures, the programs @ref build_box and @ref bin_box can
 * alternatively be used (see sections V-2.9 and V-2.11, respectively).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;topologies of single molecule for each molecule type: topo1 topo2 ...&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;coordinates of single molecule for each molecule type: pos1 pos2 ...&gt; </td></tr>
 * <tr><td> \@nsm</td><td>&lt;number of molecules for each molecule type: nsm1 nsm2 ...&gt; </td></tr>
 * <tr><td> \@dens</td><td>&lt;density of liquid (kg/m^3)&gt; </td></tr>
 * <tr><td> [\@thresh</td><td>&lt;threshold distance in overlap check; default: 0.20 nm&gt;] </td></tr>
 * <tr><td> [\@layer</td><td>&lt;create molecules in layers (along z axis)&gt;] </td></tr>
 * <tr><td> [\@boxsize</td><td>&lt;boxsize&gt;] </td></tr>
 * <tr><td> [\@fixfirst</td><td>&lt;do not rotate / shift first molecule&gt;] </td></tr>
 * <tr><td> [\@seed</td><td>&lt;random number genererator seed&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ran_box
    @topo       ic4.top   urea.top   h2o.top
    @pos        ic4.g96   urea.g96   h2o.g96
    @nsm        1         153        847
    @pdb        r     
    @dens       1000
    @thresh     0.2
    @seed       100477
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <iostream>
#include <unistd.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Physics.h"
#include "../src/gromos/Exception.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;



// declarations
bool overlap(System const & sys, double threshhold, Boundary * pbc);
int place_random(System & sys, Boundary * pbc, gsl_rng * rng, int layer = 0, int nlayer = 1);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "pos" << "nsm" << "dens" << "thresh" << "layer"
          << "boxsize" << "fixfirst" << "seed";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <topologies of single molecule for each molecule type: topo1 topo2 ...>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t@pos      <coordinates of single molecule for each molecule type: pos1 pos2 ...>\n";
  usage += "\t@nsm      <number of molecules for each molecule type: nsm1 nsm2 ...>\n";
  usage += "\t@dens     <density of liquid (kg/m^3)>\n";
  usage += "\t[@thresh   <threshold distance in overlap check; default: 0.20 nm>]\n";
  usage += "\t[@layer    (create molecules in layers (along z axis))]\n";
  usage += "\t[@boxsize  <boxsize>]\n";
  usage += "\t[@fixfirst (do not rotate / shift first molecule)]\n";
  usage += "\t[@seed     <random number genererator seed>]\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // define some variables
    //const double pi = acos(-1.0);
    const double nano_i = 1 / gmath::physConst.get_nano();
    // 1.66054 * 10^-27 * 10^27 = 1.66054
    const double fac_amu2kg = gmath::physConst.get_atomic_mass_unit() * nano_i * nano_i * nano_i ;

    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    if (args.count("seed") > 0) {
      int s = args.getValue<int>("seed", true);
      gsl_rng_set(rng, s);
    } else {
      srand(time(NULL));
    }

    // reading input and setting some values
    if (args.count("topo") != args.count("pos") || args.count("topo") != args.count("nsm")) {
      throw gromos::Exception("ran_box", "Check the number of arguments for @topo, @pos and @nsm");
    }

    if (args.count("boxsize") >= 0 && args.count("dens") >= 0) {
      throw Arguments::Exception("don't specify both boxsize and density!");
    }
    if (args.count("boxsize") == 0) {
      throw Arguments::Exception("boxsize: <length> for cubic or "
              "<len_x len_y len_z> for rectangular box!");
    }

    args.check("nsm", 1);
    vector<int> nsm;
    Arguments::const_iterator iter = args.lower_bound("nsm");
    while (iter != args.upper_bound("nsm")) {
      nsm.push_back(atoi(iter->second.c_str()));
      ++iter;
    }

    args.check("topo", 1);
    vector<string> tops;
    iter = args.lower_bound("topo");
    while (iter != args.upper_bound("topo")) {
      tops.push_back(iter->second.c_str());
      ++iter;
    }

    args.check("pos", 1);
    vector<string> insxs;
    iter = args.lower_bound("pos");
    while (iter != args.upper_bound("pos")) {
      insxs.push_back(iter->second.c_str());
      ++iter;
    }

    bool fixfirst = false;
    {
      if (args.count("fixfirst") >= 0) {
        fixfirst = true;
        if (nsm[0] != 1)
          throw Arguments::Exception("fixfirst only allowed for a single first molecule\n"
                "(just give the first system twice!)");
      }
    }

    Vec box = 0.0;
    double vtot = 0.0;
    double densit = 0.0;

    // read all topologies only to get the box length (via the mass)
    double weight = 0;
    for (unsigned int tcnt = 0; tcnt < tops.size(); tcnt++) {
      InTopology it(tops[tcnt]);
      System smol(it.system());
      for (int i = 0; i < smol.numMolecules(); i++)
        for (int j = 0; j < smol.mol(i).numAtoms(); j++)
          weight += nsm[tcnt] * smol.mol(i).topology().atom(j).mass();
    }

    if (args.count("boxsize") > 0) {

      iter = args.lower_bound("boxsize");
      {
        std::istringstream is(iter->second);
        if (!(is >> box[0]))
          throw Arguments::Exception("could not read boxsize");
      }

      ++iter;
      if (iter == args.upper_bound("boxsize")) {
        box[1] = box[2] = box[0];
      } else {
        std::istringstream is(iter->second);
        if (!(is >> box[1]))
          throw Arguments::Exception("could not read boxsize");
        ++iter;
        if (iter == args.upper_bound("boxsize"))
          throw Arguments::Exception("could not read boxsize");
        is.clear();
        is.str(iter->second);
        if (!(is >> box[2]))
          throw Arguments::Exception("could not read boxsize");
      }

      vtot = box[0] * box[1] * box[2];
      if (args["pbc"] == "t") vtot /= 2;
      densit = weight * fac_amu2kg / vtot;
    } else {
      densit = args.getValue<double>("dens", true);

      vtot = (weight * fac_amu2kg) / densit;
      if (args["pbc"] == "t") vtot *= 2;
      box[0] = pow(vtot, 1.0 / 3.0);
      box[1] = box[0];
      box[2] = box[0];
    }

    double thresh = args.getValue<double>("thresh", false, 0.20);
    thresh *= thresh;

    bool layer = false;
    if (args.count("layer") >= 0) {
      layer = true;
      std::cerr << "creating molecules in layers" << std::endl;
    }

    // printing the box size
    cerr << setw(20) << "Volume :" << vtot << endl
            << setw(20) << "Mass :" << weight * fac_amu2kg << endl
            << setw(20) << "density :" << densit << endl
            << setw(20) << "cell length :" << box[0] << " x " << box[1] << " x " << box[2] << endl
            << setw(20) << "PBC : " << args["pbc"] << endl;

    // now we do the whole thing
    // new system and parse the box sizes
    System sys;
    sys.box().K()[0] = box[0];
    sys.box().L()[1] = box[1];
    sys.box().M()[2] = box[2];

    // parse boundary conditions
    Boundary *pbc;
    if (args["pbc"] == "t")
      pbc = new TruncOct(&sys);
    else
      pbc = new RectBox(&sys);

    // loop over the number of topologies.
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
      //  with origin at cog
      if (tcnt != 0 || (!fixfirst))
        fit::PositionUtils::shiftToCog(&smol);

      //loop over the number of desired mol    
      for (unsigned int i = 0; i < unsigned(nsm[tcnt]); ++i) {
        for (int moltop = 0; moltop < smol.numMolecules(); ++moltop) {

          sys.addMolecule(smol.mol(moltop));
          // no checks, rotation on first system... (anyway?)
          if (tcnt == 0 && fixfirst) continue;

          //save position of old molecule
          Molecule oldmol = smol.mol(moltop);

          do {
            //reset molecule position
            for (int p = 0; p < oldmol.numAtoms(); p++) {
              sys.mol(sys.numMolecules() - 1).pos(p) = oldmol.pos(p);
            }

            if (layer) {
              place_random(sys, pbc, rng, tcnt, tops.size());
            } else {
              place_random(sys, pbc, rng);
            }
          } while (overlap(sys, thresh, pbc));

          cerr << (i + 1) << " of " << nsm[tcnt]
                  << " copies of molecule " << tcnt + 1
                  << " already in the box. (Total number of molecules = "
                  << sys.numMolecules() << ")." << endl;
        }
      }
      cerr << "Box now with: " << sys.numMolecules() << " molecules" << endl;
    }

    // add some zero velocitie if there is no velocity at all
    // (only if velocities should be written)
    bool numPosVel = false;
    for (int mol = 0; mol < sys.numMolecules(); mol++) {
      if (sys.mol(mol).numPos() != sys.mol(mol).numVel()) {
        sys.mol(mol).initVel();
        numPosVel = true;
      }
    }
    if (numPosVel) {
      cerr << "WARNING: Some missing velocities were set to zero.\n";
    }

    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << sys;

    gsl_rng_free(rng);
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}

/**
 * checks for overlap of all molecules os sys
 * with the last molecule of sys
 */
bool overlap(System const & sys, double threshhold, Boundary * pbc) {
  if (sys.numMolecules() == 1) return false;

  const int mol2 = sys.numMolecules() - 1;

  for (int mol1 = 0; mol1 < sys.numMolecules() - 1; ++mol1) {
    for (int a1 = 0; a1 < sys.mol(mol1).numAtoms(); ++a1) {
      for (int a2 = 0; a2 < sys.mol(mol2).numAtoms(); ++a2) {

        if ((sys.mol(mol1).pos(a1) -
                (pbc->nearestImage(sys.mol(mol1).pos(a1),
                sys.mol(mol2).pos(a2),
                sys.box()))).abs2()
                < threshhold) {
          return true;
        }
      }
    }
  }
  return false;
}

int place_random(System & sys, Boundary * pbc, gsl_rng * rng, int layer, int nlayer) {
  const int mol = sys.numMolecules() - 1;

  Vec rpos;
  const Vec box_mid = 0.5 * (sys.box().K() + sys.box().L() + sys.box().M());

  while (true) {
    double r = gsl_rng_uniform(rng);
    rpos[0] = r * sys.box().K().abs();
    r = gsl_rng_uniform(rng);
    rpos[1] = r * sys.box().L().abs();
    r = gsl_rng_uniform(rng);
    rpos[2] = layer * (sys.box().M().abs() / nlayer) + r * sys.box().M().abs() / nlayer;

    // correcting rpos for pbc (truncated octahedron / triclinic)
    if (pbc->inBox(box_mid, rpos, sys.box())) break;
  }

  /**generate random rotation
   * 1. create 2 random vectors on a sphere, v1,v2
   * 2. use v1 as the new x-axis x'
   * 3. use v2-v2*v1 as the new y-axis y'
   * 4. find z' as orthogonal to x' and y'
   * rotation matrix is then given as R = [x' y' z']
   * also make sure that det(R) = 1 (-1 corresponds to a mirroring) by changing direction of z'
   */
  const double std_dev = 1.0;
  Vec v1(gsl_ran_gaussian(rng, std_dev), gsl_ran_gaussian(rng, std_dev), gsl_ran_gaussian(rng, std_dev));
  Vec v2;

  // Now make v1 unit length
  v1 /= v1.abs();
  do {
    v2 = Vec(gsl_ran_gaussian(rng, std_dev), gsl_ran_gaussian(rng, std_dev), gsl_ran_gaussian(rng, std_dev));
    // Now get orthogonal part of v2 and make of unit length
    v2 = v2 - v1.dot(v2) * v1;
    v2 /= v2.abs();
  } while (v2.abs() < 1e-5);

  // Finally get a vector orthogonal to v1 and v2'
  // can be found by a cross product of v1 and v2'

  Vec v3;
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = -v1[0] * v2[2] + v1[2] * v2[0];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  //Finally construct rotation matrix

  Matrix m(v1, v2, v3);


  //if det(m) =-1 change direction of v3 (to ensure real rotation)
  //by construction, det(m) >0, so doesn't need to be checked

  PositionUtils::rotate(sys.mol(mol), m);
  PositionUtils::translate(sys.mol(mol), rpos);

  return 0;
}
