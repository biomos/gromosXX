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
 * @file rmsd.cc
 * calculates atom-positional root-mean-square deviations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor rmsd
 * @section rmsd calculates atom-positional root-mean-square deviations
 * @author @ref rb
 * @date 26.7.06
 *
 * The structural deformation of a molecule with respect to a reference
 * structure can be expressed in terms of a root-mean-square deviation (rmsd)
 * of the position of selected atoms. Program rmsd calculates the rmsd over a 
 * molecular trajectory. If requested a least-square rotational fit is performed before
 * to the rmsd calculation. The fit
 * can be performed using a different set of atoms than the calculation of the 
 * rmsd. If no fit is required "no" should be given.
 * Optionally a separate topology can be given for the reference. The atom 
 * selections then have to be chosen such that equivalent atoms in the system and 
 * reference system are selected and that they are in the same order. To be able to 
 * check this, option \@printatoms will print the selected atoms for each topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomsrmsd</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsd&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> [\@reftopo</td><td>&lt;molecular topology file for the reference&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates (if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> [\@printatoms</td><td>print list of selected atoms] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ditrans
    @topo       ex.top
    @reftopo    ref.top
    @pbc        r cog
    @refpbc     r cog
    @time       0 0.1
    @atomsrmsd  1:CA
    @atomsfit   1:CA,C,N
    @ref        exref.coo
    @traj       ex.tr

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/gmath/Matrix.h"
#include "../src/fit/FastRotationalFit.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Value.h"


using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;

double rmsdproperty(const utils::PropertyContainer &prop_ref, 
                             const utils::PropertyContainer &prop_sys ){

  double rmsd2=0;

  for(unsigned int i=0; i < prop_ref.size(); i++){
    
    utils::Value res = abs2(prop_ref[i]->calc() - prop_sys[i]->calc());
    rmsd2 += res.scalar();
  }
  
  return utils::sqrt(rmsd2/prop_ref.size());
}

int main(int argc, char **argv){
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atomsfit" << "atomsrmsd" << "prop" << "pbc" << "ref"
         << "time"  << "reftopo" << "refpbc" << "printatoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gathermethod>]\n";
  usage += "\t@time       <time and dt>\n";
  usage += "\t@atomsrmsd  <atoms to consider for rmsd>\n";
  usage += "\t[@atomsfit  <atoms to consider for fit>]\n";
  usage += "\t[@reftopo       <molecular topology file for the reference>\n";
  usage += "\t[@refpbc       <pbc and gather method for the reference if @reftopo is given, default: box from @ref, no gathering>\n";
  usage += "\t[@prop      <properties>\n";
  usage += "\t[@ref        <reference coordinates (if absent, the first frame of @traj is reference)>]\n";
  usage += "\t[@printatoms     < print list of selected atoms >]\n";
  usage += "\t@traj       <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);    

    System sys(it.system());
    System refSys;
    if (args.count("reftopo")>0) {
      if (args.count("ref")<=0) {
      throw gromos::Exception("rmsd",
                    "reftopo only makes sense when @ref is given");        
      }
      cerr << "# Be careful: using @reftopo will only work correctly if the atomspecifiers\n"
           << "# in atomsrmsd or atomsfit specify the equivalent atoms in \n"
           << "# both topologies in the same order!\n#\n";
      InTopology refit(args["reftopo"]);
      refSys=refit.system();
    } else {
      refSys=it.system();
    }   
     
    //print selected atoms or not
    bool printatoms=false;
    if (args.count("printatoms") >= 0) printatoms=true;
    
    // read reference coordinates...
    InG96 ic;
    if(args.count("ref")>0)
      ic.open(args["ref"]);
    else
      if(args.count("traj")>0)
	ic.open(args.lower_bound("traj")->second);

    ic.select("ALL");
    ic >> refSys;
    ic.close();
    
    Boundary *refpbc = BoundaryParser::boundary(refSys);
    Boundary::MemPtr refgathmethod=&Boundary::nogather;
    if(args.count("refpbc")>0) {
       // Parse boundary conditions and gathermethod for reference
       refpbc = BoundaryParser::boundary(refSys, args, "refpbc");
       refgathmethod= args::GatherParser::parse(refSys,refSys,args, "refpbc");
       (*refpbc.*refgathmethod)();
       // if no reference topology has been given, use normal @pbc flag
    } else if (args.count("reftopo")<=0) {
       refpbc = BoundaryParser::boundary(refSys, args);
       refgathmethod= args::GatherParser::parse(refSys,refSys,args);
       (*refpbc.*refgathmethod)();  
    } else {
      cerr << "# Warning: @reftopo but no @refpbc : reference will not be gathered\n";
    }

    // this always goes wrong. check that there is a box block in the refSys
    //if(refSys.hasBox == false && pbc->type()!='v')
    //  throw gromos::Exception("rmsd",
//			      "If pbc != v you have to give a box block "
//			      "in the reference system as well.");
    // and that we managed to read coordinates from it
    if(!refSys.hasPos)
      throw gromos::Exception("rmsd",
                              "Unable to read POSITION(RED) block from "
			      "reference positions file.");


    // System for calculation
    AtomSpecifier rmsdatoms(sys);
    AtomSpecifier rmsdatomsref(refSys);
    AtomSpecifier fitatoms(sys);
    AtomSpecifier fitatomsref(refSys);

    //get rmsd atoms
    if(args.count("atomsrmsd")>0){    
       Arguments::const_iterator iter = args.lower_bound("atomsrmsd");
       Arguments::const_iterator to = args.upper_bound("atomsrmsd");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        rmsdatoms.addSpecifier(spec);
        rmsdatomsref.addSpecifier(spec);
       }  
    }

    //try for fit atoms
    if(args.count("atomsfit") > 0){
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");

      for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        fitatomsref.addSpecifier(spec);
        fitatoms.addSpecifier(spec);
      }
    } else {
      cout << "# @atomsrmsd atoms are taken for fit." << endl;
      fitatoms=rmsdatoms;
      fitatomsref=rmsdatomsref;
    }
    
    // optionally print atomspecs to check them
    unsigned int maxsize=0;
    if (fitatoms.size() > rmsdatoms.size()) maxsize=fitatoms.size();
    else maxsize=rmsdatoms.size();
    if (printatoms) {
        for (unsigned int i=0; i < maxsize; i++) {
              cerr << "# " ;
            if (i<fitatoms.size()) {
              cerr << fitatoms.toString(i) << " "<< fitatoms.mol(i)+1 <<":"<< fitatoms.resname(i) << fitatoms.resnum(i)+1<<":" << fitatoms.name(i);
            } 
            cerr << "\t-\t";
            if (i<fitatomsref.size()) {
              cerr << fitatomsref.toString(i) << " "<< fitatomsref.mol(i)+1 <<":"<< fitatomsref.resname(i) << fitatomsref.resnum(i)+1<<":" << fitatomsref.name(i);
            }
            cerr <<  endl;
        }
    }
    
    if (fitatomsref.size() != fitatoms.size() || rmsdatomsref.size() != rmsdatoms.size())
      throw gromos::Exception("rmsd", "Number of atoms specified by atomsfit or atomsrmsd is not the same in topo and reftopo! \n@printatoms can help you find the problem.");

    // Parse boundary conditions for sys
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // the second parse argument is only used for generating the primlist (for list gathering)
    // also  in the reference; since I'm gathering the reference separately
    // I don't need it
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,sys,args);

    // Property container also for reference system
     PropertyContainer prop_ref(refSys, refpbc);
     PropertyContainer prop_sys(sys, pbc);
    {
      std::string prop;
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      // we read in all properties specified by the user
      
      for(; iter!=to; iter++){
	string spec=iter->second.c_str();
	prop += " " + spec;	
      }
      
      prop_ref.addSpecifier(prop);
      prop_sys.addSpecifier(prop);
    }


   FastRotationalFit frf;
   if(args.count("atomsrmsd") < 0 && args.count("prop") < 0){
    throw gromos::Exception("rmsd", "No rmsd atoms or property specified!");}

    int numFrames = 0;

   
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
        Matrix rot(3,3,0.0);
        for(int i=0;i<3;++i){ // initialize as identity matrix
          rot(i,i)=1.0;
        }
          
          numFrames++;
          ic.select("ALL");
	      ic >> sys >> time;
        if(!sys.hasPos)
          throw gromos::Exception("rmsd",
                             "Unable to read POSITION(RED) block from "
                              "trajectory file.");
	
	    (*pbc.*gathmethod)();

        // using FastRotationalFit::fit and rmsd methods
        // they take position vectors as input
        // shift to cog before rotational fit
        Vec cogfitref (0,0,0);
        Vec cogfitsys (0,0,0);
        if (fitatomsref.size()){
          for(unsigned int i=0;i<fitatomsref.size();++i) {
            cogfitref+=fitatomsref.pos(i);
            cogfitsys+=fitatoms.pos(i);
          }
          cogfitref/=fitatomsref.size();
          cogfitsys/=fitatomsref.size();
        }
        
        std::vector<gmath::Vec> fitref;
        std::vector<gmath::Vec> fitsys;
        std::vector<gmath::Vec> rmsdref;
        std::vector<gmath::Vec> rmsdsys;
        if (fitatomsref.size()){
          for(unsigned int i=0;i<fitatomsref.size();++i) {
            fitref.push_back(fitatomsref.pos(i)-cogfitref);
            fitsys.push_back(fitatoms.pos(i)-cogfitsys);
          }
          for(unsigned int i=0;i<rmsdatomsref.size();++i) {
            rmsdref.push_back(rmsdatomsref.pos(i)-cogfitref);
            rmsdsys.push_back(rmsdatoms.pos(i)-cogfitsys);
          }
          frf.fit(rot, fitref, fitsys);
        }

        double r = frf.rmsd(rot, rmsdref, rmsdsys);
        double rprop = rmsdproperty(prop_ref, prop_sys);

        if (args.count("atomsrmsd") > 0 && args.count("prop") > 0) {
          cout.precision(2);
          cout << time;
          cout.precision(9);
          cout << setw(15) << r  << setw(15) << rprop << endl;
        } else if(args.count("atomsrmsd") > 0){
          cout.precision(2);
          cout << time;
          cout.precision(5);
          cout << setw(10) << r  << endl;
        } else if(args.count("prop") > 0){
          cout.precision(2);
          cout << time;
          cout.precision(9);
          cout << setw(15) << rprop << endl;
        }
      }
      ic.close();
    }

    //if (rf != NULL) {
     // delete rf->getReference();
      //delete rf;
    //}
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

