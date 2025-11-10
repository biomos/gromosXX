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
 * @file bilayer_dist.cc
 * Calculates probability distribution along bilayer normal
 */

/**
 * @page programs Program Documentation
 *
 * @anchor bilayer_dist
 * @section bilayer_dist Calculates probability distribution along bilayer normal
 * @author @ref bh
 * @date 8-7-09
 *
 * Distributions along the bilayer normal are useful to characterise membrane
 * systems. This program calculates the distribution of a given set of atoms
 * with respect to the center of mass of another given set of atoms (usually,
 * this set will include all bilayer atoms to give a distribution of atoms with
 * respect to the bilayer center of mass).
 *
 * By default, all distributions are normalised to unity. However, the user
 * may be interested in density distributions. In this case, the flag
 * \@density must be included.
 *
 * The user can also obtain mass or electron densities with the help of the
 * \@mult argument. This will multiply the distribution by a given number. If the
 * user wants to calculate electron density profiles the \@mult and \@density
 * arguments must be given.
 *
 * In some cases, the bilayer is not centered in the periodic box and the center
 * of mass will not be in between the two layers but in the bulk of the solvent
 * instead. The argument \@translate can be used to circumvent this problem.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@selection</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> [\@translate</td><td>&lt;optional to translate box&gt;] </td></tr>
 * <tr><td> [\@grid</td><td>&lt;integer; default: 100&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> [\@mult</td><td>&lt;double; default: 1&gt;] </td></tr>
 * <tr><td> [\@density</td><td>&lt;optional to calculate density&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  bilayer_dist
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @atoms            1-72:a
    @selection        1-72:P
    @translate
    @density
    @mult             15
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/groTime.h"
#include "../src/gcore/Box.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;
using namespace fit;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "atoms" << "selection" << "translate"
          << "grid" << "traj" << "mult" << "density";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@atoms          <AtomSpecifier>\n";
  usage += "\t@selection      <AtomSpecifier>\n";
  usage += "\t[@translate      (only in case com is not in between layers)]\n";
  usage += "\t[@grid          <integer; default: 100>]\n";
  usage += "\t@traj           <trajectory files>\n";
  usage += "\t[@mult          <double; default: 1>]\n";
  usage += "\t[@density       (in this case the density is calculated)]";


  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    utils::Time time(args);

    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms");
      for(Arguments::const_iterator iter = args.lower_bound("atoms"); iter != to; iter++)
        atoms.addSpecifier(iter->second);
    }

    if (atoms.empty())
      throw gromos::Exception(argv[0], "@atoms: no atoms selected");

    int grid = args.getValue<int>("grid", false, 100);

    // selection
    AtomSpecifier selection(sys);
    {
      Arguments::const_iterator to = args.upper_bound("selection");
      for(Arguments::const_iterator iter = args.lower_bound("selection"); iter != to; iter++)
        selection.addSpecifier(iter->second);
    }

    // Check if translate
    bool translate = false;
    if(args.count("translate") >=0) translate = true;

    // Check if density
    bool density = false;
    if(args.count("density") >=0) density = true;

    // mult
    double mult = args.getValue<double>("mult", false, 1.0);
    
    
    vector<double> zcoord;

    double min_d;
    double max_d;
    
    // loop over all trajectories
    InG96 ic;

    int frames = 0;
    cout << "# DISTRIBUTION: " << endl;
    
    for(Arguments::const_iterator iter = args.lower_bound("traj"), to = args.upper_bound("traj");
      iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while(!ic.eof()) {
        ic.select("ALL");
        ic >> sys >> time;
        // check selection
        if (selection.empty())
          throw gromos::Exception(argv[0], "@selection: argument is empty!");

        // gather system
        (*pbc.*gathmethod)();

        Vec cm = PositionUtils::com(sys,atoms);
        

        if (translate == true){
          cm[2] = cm[2] - sys.box().M().abs()/2; 
        }
        
        min_d = -sys.box().M().abs()/2; // before was: cm[2]-sys.box()[2]/2;
        max_d = sys.box().M().abs()/2;
        
        for(unsigned int i = 0; i < selection.size(); i++) {
          Vec & vector1 = selection.pos(i);
          Vec dist = vector1 - pbc->nearestImage(vector1, cm, sys.box());
          zcoord.push_back(dist[2]);
        }
        frames++;
      } // while frames
      ic.close();

    } // for traj
    gmath::Distribution dist(min_d, max_d, grid);

    for(int unsigned i = 0; i < zcoord.size(); i++) {
      dist.add(zcoord[i]);
    }

    double binwidth = (max_d - min_d)/grid;
    
    if(density == false) {
      for(int k = 0; k < grid; k++) {
        double r = dist.value(k);
        cout << r << "\t" << double(dist[k]) * mult / (binwidth * dist.nVal()) << endl;
      }
    }
    double area = sys.box().K().abs() * sys.box().L().abs();
    if(density == true) {
      for(int k = 0; k < grid; k++) {
        double r = dist.value(k);
        cout << r << setw(15) << double(dist[k]) * mult / (area * binwidth * frames) << endl;
      }
    }
    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


