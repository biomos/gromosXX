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
 * @file iwdcf.cc
 * calculates the ion-waterdipole orientation correlation function
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor iwdcf
 * @section iwdcf calculate the ion-waterdipole orientation correlation function
 * @author @ref mk
 * @date 16. 3. 2005
 *
 * Based on the rdf program this program calculates the ion-waterdipole 
 * orientation correlation function and prints the result in form of a
 * distribution.
 * 
 * WARNING: This program works only for water, nothing else!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@centre</td><td>&lt;@ref AtomSpecifier "atoms" to take as centre&gt; </td></tr>
 * <tr><td> \@nsm</td><td>&lt;number of solvent molecules; </td></tr>
 * <tr><td> \@with</td><td>&lt;@ref AtomSpecifier "atoms" to calculate distances for&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;maximum distance&gt; </td></tr>
 * <tr><td> \@grid</td><td>&lt;number of points&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rdf
    @topo   ex.top
    @pbc    r
    @centre 1:45
    @nsm    7000
    @with   s:OW
    @cut    3.0
    @grid   100
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/Reference.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "centre" << "with" << "cut" << "grid" << "nsm"
          << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@centre <type> or\n";
  usage += "\t        <atoms> or\n";
  usage += "\t        <cog> or\n";
  usage += "\t        all\n";
  usage += "\t@with   <type> or\n";
  usage += "\t        <atoms> or\n";
  usage += "\t        all\n";
  usage += "\t@nsm    <number of solvent molecules>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t@traj   <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // read in number of solvent molecules
    int nsm = args.getValue<int>("nsm");

    // set centre atoms
    int sol_c = 0;

    AtomSpecifier centre(sys);

    {
      Arguments::const_iterator iter = args.lower_bound("centre");
      Arguments::const_iterator to = args.upper_bound("centre");
      int error = 1;

      if (iter != to) {
        string s = iter->second.c_str();
        iter++;
        if (s == "type") {
          error = 0;
          for (; iter != to; iter++) {
            string name = iter->second.c_str();
            for (int i = 0; i < sys.numMolecules(); i++)
              for (int j = 0; j < sys.mol(i).topology().numAtoms(); j++)
                if (name == sys.mol(i).topology().atom(j).name())
                  centre.addAtom(i, j);
            for (int i = 0; i < sys.sol(0).topology().numAtoms(); i++)
              if (name == sys.sol(0).topology().atom(i).name())
                for (int j = 0; j < nsm; j++) {
                  int off = j * sys.sol(0).topology().numAtoms();
                  centre.addAtom(-1, i + off);
                  sol_c++;
                }

          }

        }
        if (s == "atom") {
          error = 0;
          for (; iter != to; iter++) {
            string spec = iter->second.c_str();
            centre.addSpecifier(spec);
          }

        }
        if (s == "cog") {
          error = 0;
          AtomSpecifier temp(sys);
          centre.addAtom(-2, 0);


          for (; iter != to; iter++) {
            string spec = iter->second.c_str();
            centre.addSpecifier(spec);
          }
        }
        if (s == "all") {
          error = 0;
          for (int i = 0; i < sys.numMolecules(); i++)
            for (int j = 0; j < sys.mol(i).numAtoms(); j++)
              centre.addAtom(i, j);
        }
        if (error == 1 || centre.size() == 0)
          throw gromos::Exception("iwdcf @centre", s +
                " unknown or no atoms specified. Give 'type', 'atom', 'cog' or 'all'");
      }
    }
    // set atom to consider
    AtomSpecifier with(sys);
    int sol_w = 0;

    {
      Arguments::const_iterator iter = args.lower_bound("with");
      Arguments::const_iterator to = args.upper_bound("with");

      if (iter != to) {
        string s = iter->second.c_str();
        iter++;
        int error = 1;

        if (s == "type") {
          error = 0;

          for (; iter != to; iter++) {
            string name = iter->second.c_str();
            for (int i = 0; i < sys.numMolecules(); i++)
              for (int j = 0; j < sys.mol(i).topology().numAtoms(); j++)
                if (name == sys.mol(i).topology().atom(j).name())
                  with.addAtom(i, j);
            for (int i = 0; i < sys.sol(0).topology().numAtoms(); i++)
              if (name == sys.sol(0).topology().atom(i).name())
                for (int j = 0; j < nsm; j++) {
                  int off = j * sys.sol(0).topology().numAtoms();
                  with.addAtom(-1, i + off);
                  sol_w++;
                }
          }

        }
        if (s == "atom") {
          error = 0;

          for (; iter != to; iter++) {
            string spec = iter->second.c_str();
            with.addSpecifier(spec);
          }

        }
        if (s == "all") {
          error = 0;

          for (int i = 0; i < sys.numMolecules(); i++)
            for (int j = 0; j < sys.mol(i).numAtoms(); j++)
              with.addAtom(i, j);
        }
        if (error == 1 || with.size() == 0)
          throw gromos::Exception("iwdcf @with", s +
                " unknown or no atoms specified.\nGive 'type', 'atom' or 'all'." +
                "(is nsm=0 ?)");
      }
    }

    // read in cut-off distance
    double cut = args.getValue<double>("cut", false, 1.0);

    // read in grid number
    int grid = args.getValue<int>("grid", false, 100);

    // Parse boundary conditions
    double vol_corr = 1;

    Boundary *pbc;
    try {
      char b = args["pbc"].c_str()[0];
      switch (b) {
        case 't':
          pbc = new TruncOct(&sys);
          vol_corr = 0.5;
          break;
        case 'v':
          pbc = new Vacuum(&sys);
          break;
        case 'r':
          pbc = new RectBox(&sys);
          break;
        default:
          throw gromos::Exception("Boundary", args["pbc"] +
                  " unknown. Known boundaries are t, r and v");

      }
    } catch (Arguments::Exception &e) {
      pbc = new Vacuum(&sys);
    }


    for (int i = 0; i < with.size(); i++) {
      if (with.mol(i) >= 0)
        throw gromos::Exception("iwdcf",
              " works only for water, i.e. solvent...Abort!");
    }


    // define input coordinate
    InG96 ic;

    // set up distribution arrays

    std::vector<double> iwdcf(grid);
    gmath::Distribution dist(0, cut, grid);


    for (int i = 0; i < grid; i++) iwdcf[i] = 0;

    // loop over all trajectories
    int count_frame = 0;

    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      if (sol_c || sol_w) ic.select("ALL");


      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys;

        if (nsm > sys.sol(0).numPos() / sys.sol(0).topology().numAtoms())
          throw gromos::Exception("iwdcf",
                " nsm specified is more than in coordinate file");
        else
          sys.sol(0).setNumPos(nsm * sys.sol(0).topology().numAtoms());

        //pbc->gather();
        // loop over the centre atoms
        int start = 0;

        Vec cog(0.0, 0.0, 0.0);

        if (centre.mol(0) == -2) {
          start = centre.size() - 1;
          for (int i = 1; i < centre.size(); i++)
            cog += sys.mol(centre.mol(i)).pos(centre.atom(i));
          cog /= (centre.size() - 1);
        }


        // now really loop over the centre atoms
        for (int i = start; i < centre.size(); i++) {

          //gmath::Distribution distM(0, cut, grid);

          Vec curr;
          if (centre.mol(i) >= 0)
            curr = sys.mol(centre.mol(i)).pos(centre.atom(i));
          else
            curr = sys.sol(0).pos(centre.atom(i));
          // see if this atom is also in the with-list
          int inwith = 0;

          for (int j = 0; j < with.size(); j++)
            if (with.atom(j) == centre.atom(i) && with.mol(j) == centre.atom(j))
              inwith = 1;

          if (centre.mol(0) == -2) curr = cog;

          //loop over the atoms to consider
          for (int j = 0; j < with.size(); j++) {
            //calculate distance only if this atom is not the current centre
            if (!(with.mol(j) == centre.mol(i) && with.atom(j) == centre.atom(i))) {
              Vec tmp;

              tmp = pbc->nearestImage(curr,
                      sys.sol(0).pos(with.atom(j)),
                      sys.box());
              Vec Ionwater(tmp - curr);

              double dis = Ionwater.abs();
              dist.add(dis);

              Vec H1;
              Vec H2;

              H1 = pbc->nearestImage(tmp,
                      sys.sol(0).pos(with.atom(j) + 1),
                      sys.box());

              H2 = pbc->nearestImage(tmp,
                      sys.sol(0).pos(with.atom(j) + 2),
                      sys.box());


              Vec t1(H1 - tmp);
              Vec t2(H2 - tmp);



              Vec e(t1 + t2);
              e /= e.abs();

              if (dist.inrange(dis)) iwdcf[dist.getbin(dis)] += ((Ionwater.dot(e)) / dis);


            }
          }

        }
        count_frame++;
      }
      ic.close();
    } //end loop over trajectory
    //now correct the distribution for the number of frames and the number 
    //of centre atoms
    cout << "# number of frames considered: " << count_frame << endl;
    int divide = count_frame;

    if (centre.mol(0) != -2) divide *= centre.size();

    for (int i = 0; i < grid; i++) {
      double r = (double(i) + 0.5) * cut / grid;
      double out = 0;
      if (dist[i] == 0 && iwdcf[i] == 0) out = 0;
      else out = iwdcf[i] / ((double(divide) * double(dist[i])));
      cout << r << "\t" << out << endl;
    }

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

