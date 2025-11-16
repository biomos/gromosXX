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
 * @file structal.cc
 * Tool to compare two structures
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor structal
 * @section structal threads a topology through a structure
 * @author @ref as @ref co
 * @date 20-12-10
 *
 * This program can be used to compare the structure of two molecules regardless
 * of their sequence.
 * You will need:
 *      - a molecular topology of the backbone of your reference structure
 *      - a molecular topology and coordinates of the backbone of the structure you want the
 *        reference structure to be compared to
 *
 * The program will thread the reference topology over the topology to compare to
 * and write out a trajectory containing the coordinates from the system to
 * compare to but matching the reference topology.
 * This trajectory can then be analysed with the program @ref rmsd using the 
 * orginal coordinates from the reference structure as reference.
 * 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@reftopo</td><td>&lt; reference molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@select</td><td>&lt;@ref AtomSpecifier "residues" to consider&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> [\@verbose]</td><td></td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 structal
    @topo      ex.top
    @reftopo   refex.top
    @pbc       r
    @pos       ex.g96
    @select    1:res(1-5:a)
    @verbose
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/OutformatParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gcore/Box.h"
#include "../src/gromos/Exception.h"



using namespace args;
using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;



int main(int argc, char **argv){

// Argument list
  Argument_List knowns;
  knowns << "topo" << "reftopo" << "pbc" << "pos" << "select" << "verbose";

// usage
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file that you want to compare to>\n";
  usage += "\t@reftopo     <molecular topology file that you want to compare>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@pos         <coordinate file that you want to compare to>\n";
  usage += "\t@select      <residues to compare (e.g. 1:res(1-5:a)>\n";
  usage += "\t[@verbose]\n";

// try to do something
  try{
// start code //
    Arguments args(argc, argv, knowns, usage);

    // debug
    bool debug = false;
    if(args.count("verbose")>=0) {
        debug = true;
    }
    // end debug

     if (debug) {
         // start
         cerr << "Starting..." << endl;
     }

    // get topo 1
     InTopology it1(args["reftopo"]);
     System sys1(it1.system());

     if (debug) {
         cerr << "Reference topology read in (sys1)..." << endl;
         cerr << "\tContains " << sys1.numMolecules() << " molecule";
                 if(sys1.numMolecules() > 1) cerr << "s";
         cerr << endl;
     }

     // prepare coordinates for sys1
     for(int m=0; m < sys1.numMolecules(); m++)
         sys1.mol(m).initPos();
     if (debug) {
         cerr << "Initialized coordinates for sys1..." << endl;
     }
     
     // get topo 2
     InTopology it2(args["topo"]);
     System sys2(it2.system());
     System refSys2(it2.system());

     if (debug) {
         cerr << "Topology read in (sys2)..." << endl;
         cerr << "\tContains " << sys2.numMolecules() << " molecule";
                 if(sys2.numMolecules() > 1) cerr << "s";
         cerr << endl;
     }
 
     // parse boundary conditions for the second system
        Boundary *pbc = BoundaryParser::boundary(sys2, args);
        Boundary::MemPtr gathmethod = args::GatherParser::parse(sys2,refSys2,args);

     if (debug) {
         cerr << "Parsed pbc" << endl;
     }

      // read in coords into second system
      InG96 ic(args["pos"]);
      ic.select("SOLUTE");
      ic >> sys2;
      ic.close();
      // apply pbc
     (*pbc.*gathmethod)();

      if (debug) {
         cerr << "Read in coordinates (sys2) from " <<  args["pos"] << endl;
      }

     // create specifierers
      utils::AtomSpecifier ls1(sys1);
      {
      Arguments::const_iterator iter=args.lower_bound("select"),
	to=args.upper_bound("select");
      for(;iter!=to; ++iter)
	ls1.addSpecifier(iter->second);
      }

      // get the total atom count including gaps
      // we start from the first atom regardless of ls1, because in the copying
      // that is the reference (atomcount + j + lcount)
      ls1.sort();
      //int firstm=ls1.mol(0);
      //int firsta=ls1.atom(0);
      int firstm=0;
      int firsta=0;
      int lastm=ls1.mol(ls1.size()-1);
      int lasta=ls1.atom(ls1.size()-1);
      int numawg=0;
      for(int m=firstm; m <= lastm; m++) {
          numawg+=sys1.mol(m).numAtoms();
      }
      numawg -= (sys1.mol(lastm).numAtoms() - lasta - 1 + firsta);

      utils::AtomSpecifier ls2(sys2);
      ls2.addSpecifier("a:a");

      if (debug) {
        cerr << "Created Atomspecifiers (sys1 and sys2)..." << endl;
         cerr << "\tsys1: " << args["reftopo"] << endl;
         cerr << "\t\t Contains: " << ls1.mol(ls1.size()-1)+1 << " molecule";
                 if(ls1.mol(ls1.size()-1)+1 > 1) cerr << "s";
         cerr << endl;
         cerr << "\t\t Contains: " << ls1.size() << " atoms"<< endl;
         cerr << "\t\t Contains: " << numawg << " atoms with gaps"<< endl;
         cerr << "\tsys2: " << args["topo"] << endl;
         cerr << "\t\t Contains: " << ls2.mol(ls2.size()-1)+1 << " molecule";
                 if(ls2.mol(ls2.size()-1)+1 > 1) cerr << "s";
         cerr << endl;
         cerr << "\t\t Contains: " << ls2.atom(ls2.size()-1)+1 << " atoms"<< endl;
      }

      if (debug) {
       cerr << "Looping with reftopo over topo..." << endl;
      }

      //
      // output
      //
      // define and open
      OutG96 oc;
      oc.open(cout);
      // write title
      std::ostringstream title;
      title << "Stepwise structural comparison" << endl;
      title << "Comparing "  << args["reftopo"] << " to " << args["topo"] << endl;
      oc.writeTitle(title.str());
  //
  // diffrent format
  //   string ext;
  //  OutCoordinates & oc = *OutformatParser::parse(args, ext);
  //
      // end output defination

      //
      // loop
      //
         //
         // initate counters
         //
            // how many atoms do we have to compare to
            int sysna = ls2.size();
            // how many atoms has the reference system
            // numawg
            // whats the end
            int lend = sysna - numawg;
            // loop counter
            int lcount = 0;
            // timestep
            int ltime = 0;
            double ldouble = 0;

         // residue numbers
         //   int sysnr = 0;
         //   for(int m=0; m < sys2.numMolecules(); m++){
         //      sysnr += sys2.mol(m).topology().numRes();
         //  }
         //  int rcount=0;

      // begin loop
      for(;lcount<=lend; ){

      // copy coordinates
      for(int i=0, j=0; i<ls1.size(); i++) {
          if(ls1.name(i) == ls2.name(ls1.atom(i)+j+lcount)) {
             ls1.pos(i) = ls2.pos(ls1.atom(i)+j+lcount);
          } else {
              j++;
              i--;
         }
      }
      // copy box
      sys1.box() = sys2.box();

      // write output
      // title
      oc.writeTimestep(ltime, ldouble);
      oc << ls1;

      // update counters
        // put time one forward
        ltime = ltime + 1;
        ldouble = ldouble + 1;
        // go one residue forward
        int a=0;
        while(ls2.resnum(lcount)==ls2.resnum(lcount+a)) {
            a++;
        }
        lcount += a;
        if (debug) {
        cerr << "This is step " << ltime << "\r";
        }

     }
     //end loop
     if (debug) {
     // end the counter line
     cerr << endl;
     // done
     cerr << "Done" << endl;
     }

     // close outstream
     oc.close();

  // throw something
  // throw gromos::Exception("thrown", "outch");
   
// end code //
  }
// catch all
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
