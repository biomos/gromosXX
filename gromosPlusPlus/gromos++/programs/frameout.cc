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
 * @file frameout.cc
 * Write out individual snapshots of a molecular trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor frameout
 * @section frameout Write out individual snapshots of a molecular trajectory
 * @author @ref mk
 * @date 24.8.06
 *
 *
 * Program frameout can be used to extract individual snapshots from a
 * molecular trajectory file. Different formats are supported e.g. the
 * GROMOS format, the PDB format and an VMD-Amber format which can be read
 * in by VMD. The user determines which frames should be written out
 * and if solvent should be included or not. Atom positions can be corrected
 * for periodicity by taking the nearest image to connected atoms, or to the
 * corresponding atom in a reference structure. A least-squares rotational fit
 * to a reference structure can be performed based on selected atoms.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> [\@spec</td><td>&lt;specification for writing out frames: ALL (default), EVERY or SPEC&gt;] </td></tr>
 * <tr><td> [\@frames</td><td>&lt;frames to be written out&gt;] </td></tr>
 * <tr><td> \@outformat</td><td>&lt;@ref args::OutformatParser "output coordinates format"&gt; </td></tr>
 * <tr><td> [\@include</td><td>&lt;SOLUTE (default), SOLVENT or ALL&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference structure to fit to&gt;] </td></tr>
 * <tr><td> [\@gathref</td><td>&lt;reference structure to gather with respect to(use ggr as gather method)&gt;] </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to fit to&gt;] </td></tr>
 * <tr><td> [\@name</td><td>&lt;prefix for the filenames&gt;] </td></tr>
 * <tr><td> [\@single</td><td>&lt;write to a single file&gt;] </td></tr>
 * <tr><td> [\@notimeblock</td><td>&lt;do not write timestep block&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  frameout
    @topo        ex.top
    @pbc         r
    @spec        SPEC
    @frames      1 3
    @outformat   pdb
    @include     ALL
    @ref         exref.coo
    @gathref     exref.coo
    @atomsfit    1:CA
    @name        protein
    [@notimeblock   ]
    [@time        0 2]
    @single
    @traj        ex.tr
 @endverbatim
 *
 * <hr>
 */
// frameout.cc

#include <cassert>
#include <cctype>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

bool writeFrame(int i, std::string const & spec, vector<int> const & fnum,
        unsigned int & framesWritten, bool & done);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "pbc" << "spec" << "frames" << "outformat"
          << "include" << "ref" << "atomsfit" << "single" << "notimeblock" << "time" << "addvirtual" << "name";
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t[@spec      <specification for writing out frames: ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@frames    <frames to be written out>]\n";
  usage += "\t[@outformat <output coordinates format>]\n";
  usage += "\t[@include   <SOLUTE (default), SOLVENT, VIRTUAL or ALL>]\n";
  usage += "\t[@ref       <reference structure to fit to>]\n";
  usage += "\t[@atomsfit  <atoms to fit to>]\n";
  usage += "\t[@name       <prefix and postfix of output trajectories>]\n";
  usage += "\t[@single    <write to a single file>]\n";
  usage += "\t[@addvirtual  <write coordinates for virtual atoms>]\n";
  usage += "\t[@notimeblock <do not write timestep block>]\n";
  usage += "\t[@time      <time and dt>]\n";
  usage += "\t@traj       <trajectory files>\n";
  //  usage += "\t[@gathref   <reference structure to gather with respect to"
  //    "(use ggr as gather method)>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // create Time object
    utils::Time time(args);

    // get @notimeblock argument
    bool notimeblock = false;
    if (args.count("notimeblock") >= 0)
      notimeblock = true;

    // get simulation time either from the user or from the files
    //bool usertime=false;
    
    //if (args.count("time") > 0) {
    //  usertime=true;
    //}

    // do we want to fit to a reference structure
    bool fit = false;
    //System refSys(sys);
    System refSys(it.system());
    Reference reffit(&refSys);
    Vec cog(0.0, 0.0, 0.0);

    // now we always define a reference
    if (args.count("ref") > 0) {
      fit = true;

      // read reference coordinates...
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> refSys;
      ic.close();
    } else {
      InG96 ic;
      if (args.count("traj") > 0) {
        ic.open(args.lower_bound("traj")->second);

        ic.select("ALL");
        ic >> refSys;
        ic.close();
      } else {
        throw gromos::Exception("frameout", "no trajectory specified (@traj)");
      }
    }

    AtomSpecifier fitatoms(refSys);

    bool name = false;
    if (args.count("name") >= 0)
      name = true;

    //try for fit atoms
    if (args.count("ref") > 0) {
      if (args.count("atomsfit") > 0) {
        Arguments::const_iterator iter = args.lower_bound("atomsfit");
        Arguments::const_iterator to = args.upper_bound("atomsfit");
        for (; iter != to; iter++) fitatoms.addSpecifier(iter->second);
      } else {
        throw gromos::Exception("frameout",
                "If you want to fit (@ref) then give "
                "atoms to fit to (@atomsfit)");
      }
      reffit.addAtomSpecifier(fitatoms);
      cog = PositionUtils::cog(refSys, reffit);
    }
    // does this work if nothing is set?
    RotationalFit rf(&reffit);

    // parse includes
    string inc = "SOLUTE";
    if (args.count("include") > 0) {
      inc = args["include"];
      transform(inc.begin(), inc.end(), inc.begin(), static_cast<int (*)(int)> (std::toupper));
      if (inc != "SOLUTE"  && inc != "ALL"  && inc != "SOLVENT" && inc != "VIRTUAL") 
        throw gromos::Exception("frameout",
              "include format " + inc + " unknown. \n"
              "Give SOLUTE, SOLVENT or ALL");
    }
    string inco = inc;
    if (args.count("addvirtual") >= 0){
      inco = inc+"V";
      // little hack: if you do want the virtual atoms, but not the solute, we
      //              still need to read in the solute
      if(inc == "SOLVENT") inc = "ALL";
    }
    if(inc == "VIRTUAL") { 
      inc = "SOLUTE";
      inco = "SOLVENTV"; 
    }

    // parse spec
    string spec = "ALL";
    vector<int> fnum;
    if (args.count("spec") > 0) {
      spec = args["spec"];
      transform(spec.begin(), spec.end(), spec.begin(), static_cast<int (*)(int)> (std::toupper));
      if (spec != "ALL" && spec != "EVERY" && spec != "SPEC")
        throw gromos::Exception("frameout",
              "spec format " + spec + " unknown. \n"
              "Give ALL, EVERY or SPEC");
      if (spec == "EVERY" || spec == "SPEC") {
        //smack in the framenumbers
        for (Arguments::const_iterator it = args.lower_bound("frames");
                it != args.upper_bound("frames"); ++it) {
          int bla = atoi(it->second.c_str());
          fnum.push_back(bla);
        }
        if (fnum.size() == 0) {
          throw gromos::Exception("frameout",
                  "if you give EVERY or SPEC you have to use "
                  "@frames as well");
        }
        if (fnum.size() != 1 && spec == "EVERY") {
          throw gromos::Exception("frameout",
                  "if you give EVERY you have to give exactly"
                  " one number with @frames");
        }
      }
    }

    // parse outformat
    bool single_file = false;
    string ext;
    OutCoordinates *oc = OutformatParser::parse(args, ext);

    //Compute the atomic radii to write them out (*.pqr format)
    utils::compute_atomic_radii_vdw(sys, it.forceField());  

    // check if single_file is overwritten by user
    if (args.count("single") >= 0 || dynamic_cast<Outvmdam*> (oc) != NULL)
      single_file = true;

    // loop over all trajectories
    InG96 ic;
    
    int numFrames = 0;
    ofstream os;

    // is output file open
    bool alopen = false;

    // number of frames that have been written.
    unsigned int framesWritten = 0;
    // all the frames are written: stop reading the topology.
    bool done = false;

    // make sure there is an argument @pbc
    args.check("pbc", 1);
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      // loop over all frames

      while (!ic.eof()) {
        numFrames++;
        ic.select(inc);
        if (!notimeblock) {
          ic >> sys >> time;
        } else {
          ic >> sys;
        }

        //pbc->setReferencefull(refSys);
        if (writeFrame(numFrames, spec, fnum, framesWritten, done)) {
          (*pbc.*gathmethod)();

          if (fit) {
            rf.fit(&sys);
            PositionUtils::translate(&sys, cog);
          }

          if ((!alopen) || (!single_file)) {
            ostringstream pdbName;
            if (name==true) {
                pdbName << std::string(args["name"])  << "_" << setw(5)<< setfill('0') << numFrames << ext;
            }
            else {
                pdbName << "FRAME_"<< setw(5)<< setfill('0') << numFrames << ext;
            }
            string file=pdbName.str();
            os.open(file.c_str());
            oc->open(os);
            oc->select(inco);    
            oc->writeTitle(file);
            alopen = true;
          }

          if (!notimeblock) { 
            oc->writeTimestep(time.steps(), time.time());
          }
          *oc << sys;

          if (!single_file)
            os.close();
        }
        if (done)
          break;
      }
      ic.close();
      if (done)
        break;
    }
    if (single_file)
      os.close();
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool writeFrame(int i, std::string const & spec, vector<int> const & fnum,
        unsigned int & framesWritten, bool & done) {
  if (spec == "ALL") {
    ++framesWritten;
    return true;
  } else if (spec == "EVERY" && i % fnum[0] == 0) {
    ++framesWritten;
    return true;
  } else if (spec == "SPEC") {
    for (unsigned int j = 0; j < fnum.size(); ++j) {
      if (fnum[j] == i) {
        ++framesWritten;
        if (framesWritten == fnum.size())
          done = true;
        return true;
      } // write frame?
    } // frames
  }
  return false;
}


