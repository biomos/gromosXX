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
 * @file bilayer_oparam.cc
 * calculates order parameters in bilayer
 */

/**
 * @page programs Program Documentation
 *
 * @anchor bilayer_oparam
 * @section bilayer_oparam calculates order parameter for lipids
 * @author @ref bh wh
 * @date 28-04-09
 *
 * Deuterium order parameters (@f$S_{CD}@f$) can be derived from deuterium
 * quadrupole splitting experiments and have being widely used to study
 * biological membranes and to validate biomolecular force- fields.
 * The corresponding carbon-hydrogen order parameters (@f$S_{CH}@f$) can be
 * calculated by computing the correlation functions describing
 * the reorientation of the carbon-hydrogen vectors. More precisely, for each
 * methylene group along the chain, an order parameter tensor
 * @f${\bf\underline {S}}@f$ can be defined as
 *
 * @f[ S_{ij}=\frac{1}{2}\langle3 cos\theta_{i}~cos\theta_{j} - \delta_{ij}\rangle @f],
 *
 * where @f$\theta_{i}@f$ is the angle between the @f$i^{th}@f$ local molecular axis
 * (@f$x'@f$, @f$y'@f$ or @f$z'@f$) and the bilayer normal (@f$z@f$-axis),
 * @f$\delta_{ij}@f$ is the
 * Kronecker delta symbol and @f$\langle@f$...@f$\rangle@f$ stands for trajectory
 * averaging. As a convention,
 * for the @f$n^{th}@f$ methylene group {\em C@f$_n@f$}, the direction of the vector
 * @f$C_{n-1}-C_{n+1}@f$ is taken as @f$z'@f$, the direction of the vector normal to
 * @f$z'@f$ in the plane @f$C_{n-1}@f$, @f$C_{n}@f$, and @f$C_{n+1}@f$ defines @f$y'@f$, while @f$x'@f$
 * is the direction of the vector perpendicular both to @f$z'@f$ and @f$y'@f$.
 * The quantity @f$S_{CH}=-(2/3S_{xx}+1/3S_{yy})@f$ is the value to be compared with
 * the experimental S@f$_{CD}@f$ value.
 *
 * The order parameters are calculated with respect to a fixed orientational vector (corresponding
 * to the direction of the experimental magnetic field; usually taken along the bilayer normal) 
 * given by the flag \@refvec.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@refvec</td><td>&lt;x y z&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  stacking
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @atoms            1-72:11-24
    @refvec           0 0 1
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

/* REMARK: This code was adapted from a pre-existing code written by
 an unknown author. The main operations are basically from the old version, but
 the code was formatted to be compatible with gromos++, documentation was added,
 some new options and variables were introduced and some clean-up was made.*/

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" 
          << "refvec" << "time"
          << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@atoms          <atomspecifier>\n";
  usage += "\t@refvec        <x y z>]\n";
  usage += "\t@traj           <trajectory files>\n";
  

  try {
    Arguments args(argc, argv, knowns, usage);

    //get reference vector, normalize it
    Vec refvec(0.0, 0.0, 0.0);

    {
      Arguments::const_iterator iter = args.lower_bound("refvec");
      if(iter != args.upper_bound("refvec")) {
        refvec[0] = atof(iter->second.c_str());
        ++iter;
      }
      if(iter != args.upper_bound("refvec")) {
        refvec[1] = atof(iter->second.c_str());
        ++iter;
      }
      if(iter != args.upper_bound("refvec")) {
        refvec[2] = atof(iter->second.c_str());
      }
    }
    //normalize
    refvec = refvec.normalize();
    // read topology
    InTopology it(args["topo"]);
    //make system out of topology
    System sys(it.system());

    // get atom specifiers
    AtomSpecifier bilayer_atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms");
      for(Arguments::const_iterator iter = args.lower_bound("atoms"); iter != to; iter++)
        bilayer_atoms.addSpecifier(iter->second);
    }


    // calculate number of molecules and number of atoms per molecule
    // NOTE: assumes same atoms specified for each molecule
    int moln = bilayer_atoms.mol(bilayer_atoms.size() - 1) + 1;
    // JRA: if you add 1, num_atperlip is wrong and you get a segmentation fault
    //int num_atperlip = (bilayer_atoms.size() + 1) / moln;
    int num_atperlip = (bilayer_atoms.size()) / moln;
    // check number of molecules
    if(moln <= 0) {
      throw gromos::Exception("bilayer_oparam", "molecule number cannot be <= 0!\n");
    }

    // store atom numbers (from first lipid) and check atoms
    vector<int> atoms;
    for(int i = 0; i < num_atperlip; i++) {
      if (bilayer_atoms.atom(i) > 0 && bilayer_atoms.atom(i) < sys.mol(moln-1).numAtoms()-1) {
      atoms.push_back(bilayer_atoms.atom(i));
      } else if (bilayer_atoms.atom(i) == 0) {
          throw gromos::Exception("bilayer_oparam", "cannot calculate the oparam for the first atom!\n");
      } else {
          throw gromos::Exception("bilayer_oparam", "cannot calculate the oparam for the last atom!\n");
      }  
    }
    // check we have some atoms
    if(int (atoms.size()) == 0) {
      throw gromos::Exception("bilayer_oparam", "at least one atom needs to be defined!\n");
    }

    // define vector
    vector<int> at;

	// get the bonded neighbors
	for(int i = 0; i < int (atoms.size()); ++i) {
	    Neighbours neighbours(sys, bilayer_atoms.mol(0), bilayer_atoms.atom(i));
		Neighbours::const_iterator itn = neighbours.begin(), ton = neighbours.end();
    // store sets of 3 atoms (bonded)
	    at.push_back(*itn);
		at.push_back(atoms[i]);
		for (;itn != ton-1; ++itn){
		} 
		at.push_back(*itn);
	}
	    

    System refSys(it.system());

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    // loop over all trajectories
    int numFrames = 0;
    Vec z(0.0, 0.0, 0.0);
    Vec y(0.0, 0.0, 0.0);
    Vec x(0.0, 0.0, 0.0);
    Vec cos(0.0, 0.0, 0.0);
    Vec scos2(0.0, 0.0, 0.0);
    vector<Vec> S;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      S.push_back(z);
    }
    vector<Vec> avcos2;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      avcos2.push_back(z);
    }
    // define input coordinate
    InG96 ic;
    for(Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()) {
        numFrames++;
        ic >> sys;
        (*pbc.*gathmethod)();

        // calculate the z-vector between atom i-1 and i+1, normalize
        //int cc = -2; // unnecessary
        for(int i = 0; i < moln; ++i) {
          int cc = -2;
          for(int j = 0; j< int (at.size() / 3); j++) {
            cc += 2;
            z = sys.mol(i).pos(at[j + 2 + cc]) - sys.mol(i).pos(at[j + cc]);
            z = z.normalize();
            //calculate y, normalize
            y = (sys.mol(i).pos(at[j + 1 + cc]) - sys.mol(i).pos(at[j + cc]));
            y = y - ((sys.mol(i).pos(at[j + 1 + cc])-(sys.mol(i).pos(at[j + cc]))).dot(z)) * z;
            y = y.normalize();
            //calculate x
            x = z.cross(y);
            x = x.normalize(); //is this necessary?
            // determine the angle
            cos[0] = refvec.dot(x);
            cos[1] = refvec.dot(y);
            cos[2] = refvec.dot(z);
            scos2 = cos*cos;
            avcos2[j] += scos2;
          }
        }

      }
    }
    ic.close();

    // average the avcos2 finally
    for(int i = 0; i < int (at.size() / 3); ++i) {
      avcos2[i] /= (numFrames * moln);
    }

    //get the S components
    for(int i = 0; i< int (at.size() / 3); ++i) {
      Vec tmp = avcos2[i];
      tmp[0] = ((3 * tmp[0]) - 1) / 2;
      tmp[1] = ((3 * tmp[1]) - 1) / 2;
      tmp[2] = ((3 * tmp[2]) - 1) / 2;
      S[i] = tmp;
    }

    //determine the  SCDOP and SCDAP, SCDEP
    vector<double> SCDOP, SCDAP, SCDEB;
    for(int i = 0; i< int (at.size() / 3); ++i) {
      Vec tmp = avcos2[i];
      Vec tmp1 = S[i];
      SCDOP.push_back(-(tmp[0]+(tmp[1] - 1) / 2));
      SCDAP.push_back(tmp1[2] / 2);
      SCDEB.push_back((2 * tmp1[0] / 3) + tmp1[1] / 3);
    }

    //print out results
    //REMARK: The old output was changed to agree with doxygen:
    //SCD -> SCH (just a matter of convention)
    //The variable names however are still the same.
    cout << setw(4) << "Atom"
            << setw(8) << "SX" << setw(8) << "SY" << setw(8) << "SZ"
            << setw(8) << "SCHOP" << setw(8) << "|SCD|"<<endl;
    // In the old version this 2 quantities were also printed to the output
    // But I don't know the utility of that... so I removed them
            //<< setw(8) << "SCHAP"
            //<< setw(8) << "SCHEB" << endl;
    cout << endl;

    int c = -1;
    for(int i = 0; i < int (at.size() / 3); ++i) {
      c += 2;
      Vec tmp = S[i];
      cout.setf(ios::floatfield, ios::fixed);
      cout.setf(ios::right, ios::adjustfield);
      cout.precision(4);
      cout << setw(4) << at[i + c] + 1
              << setw(8) << tmp[0] << setw(8) << tmp[1] << setw(8) << tmp[2]
              << setw(8) << SCDOP[i] << setw(8) << fabs(tmp[1]) << setw(8) << SCDEB[i] << endl;
    }
  }  catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


