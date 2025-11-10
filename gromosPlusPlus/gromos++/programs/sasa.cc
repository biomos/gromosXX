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
 * @file sasa.cc
 * Calculates solvent-accessible surface areas for selected atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sasa
 * @section sasa Calculates solvent-accessible surface areas for selected atoms
 * @author @ref mk
 * @date 21-6-07
 *
 * Program sasa calculates and prints the solvent-accessible surface
 * area (sasa) of selected heavy atoms (sasaatoms) in the solute part of the molecular system.
 * It also calculates the contribution made by a specified set of heavy atoms.
 * The program uses the algorithm of Lee and Richards [J. Mol. Biol., 55, 379-400 (1971)].
 * A spherical probe of given radius is rolled over the surface of the sasaatoms
 * (the size of the probe is typically 0.14 nm for water). The path traced out
 * by its centre gives the accessible surface.
 *
 * In GROMOS, the radii of the heavy atoms are obtained by calculating
 * the minimum energy distance of the interaction between the heavy
 * atom and the first solvent atom. This value is reduced by the specified
 * probe radius to account for the radius of the solvent atom.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to monitor&gt; </td></tr>
 * <tr><td> \@sasaatoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider for sasa (default: all)&gt; </td></tr>
 * <tr><td> [\@zslice</td><td>&lt;distance between the Z-slices through the molecule (default: 0.005~nm)&gt;] </td></tr>
 * <tr><td> \@probe</td><td>&lt;probe IAC and radius&gt; </td></tr>
 * <tr><td> [\@verbose</td><td>(print summaries)] </td></tr>
 * <tr><td> [\@skip</td><td>skip first n frames] </td></tr>
 * <tr><td> [\@skip</td>use only every nth frame <td>] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  sasa
    @topo     ex.top
    @pbc      r
    @time     0 1
    @atoms    1:CB
    @zslice   0.005
    @probe    5 0.14
    @verbose
    @traj     ex.tr
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
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/groTime.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace utils;

void heapsort(double* values, int n, int* key);

//some constants

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "zslice" << "atoms" << "sasaatoms" << "probe" << "traj"
          << "verbose" << "skip" << "stride";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type>\n";
  usage += "\t[@time      <time and dt>]\n";
  usage += "\t@atoms      <atoms to monitor>\n";
  usage += "\t[@sasaatoms <atoms to consider for sasa (default: all)>]\n";
  usage += "\t[@zslice    <distance between the Z-slices (default: 0.005)>]\n";
  usage += "\t@probe      <probe IAC and radius>\n";
  usage += "\t[@verbose   (print summaries)\n";
  usage += "\t[@skip      <n> skip first n frames>]\n";
  usage += "\t[@stride    <n> use only every nth frame]\n";
  usage += "\t@traj       <trajectory files>\n";

  try {

    // set some constants
    double const PI = gmath::physConst.get_pi();
    double const twoPI = 2 * PI;
    
    Arguments args(argc, argv, knowns, usage);

    //   get simulation time
    Time time(args);

    //try for zslice, else set default
    double zslice = args.getValue<double>("zslice", false, 0.005);

    //@skip and @stride
    int stride = args.getValue<int>("stride", false, 1);
    int skip = args.getValue<int>("skip", false, 0);

    //get probe IAC and radius
    int probe_iac;
    double probe;
    if (args.count("probe") > 1) {
      vector<double> probearg = args.getValues<double>("probe", 2);
      probe_iac = int(probearg[0]);
      probe = probearg[1];
    } else {
      throw gromos::Exception("sasa", "You need to specify the probe IAC and radius!");
    }

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());


    //get sasa atoms: these are the atoms to roll the probe over
    AtomSpecifier rollatoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("sasaatoms");
      Arguments::const_iterator to = args.upper_bound("sasaatoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        rollatoms.addSpecifier(spec);
      }
      if(rollatoms.size()==0){
	for(int i=0; i < sys.numMolecules(); i++){
	  for(int j=0; j < sys.mol(i).numAtoms(); j++){
	    rollatoms.addAtom(i,j);
	  }
	}
      }
    }
      
    //get atoms: these are the atoms to report separately
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atoms.addSpecifier(spec);
      }
      // check that they are all in 'roll atoms'
      for(unsigned int i=0; i < atoms.size(); ++i){
	if(rollatoms.findAtom(atoms.mol(i), atoms.atom(i))==-1){
	  ostringstream os;
	  os << "Error: not all selected atoms are part of the sasaatoms: " << atoms.toString(i);
	  throw (gromos::Exception("sasa", os.str()));
	}
      }
    }

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    //get radii and determine heavy atoms from rollatoms
    AtomSpecifier heavyatoms(sys);
    vector<double> radheavy;
    vector<double> radheavysq;
    double rmax = 0;

    utils::compute_atomic_radii_vdw(probe_iac, probe, sys, it.forceField());

    //get all heavy atoms...
    for(unsigned int i = 0; i < rollatoms.size(); ++i) {
      if (!sys.mol(rollatoms.mol(i)).topology().atom(rollatoms.atom(i)).isH()) {
	heavyatoms.addAtom(rollatoms.mol(i), rollatoms.atom(i));
	double rad = rollatoms.radius(i);
	rad += probe;
	rmax = ((rad > rmax) ? (rad) : (rmax));
	radheavy.push_back(rad);
	radheavysq.push_back(rad * rad);
      }
    }

    // define input coordinate
    InG96 ic;
    
    // declare some variables
    vector<double> accs(heavyatoms.size(), 0.0); // stores the accessibility

    // print title
    cout << "#               "
	 << setw(15) << "selected" << ' '
	 << setw(15) << "heavy sasa" << endl
	 << "# time          "
	 << setw(15) << "atoms" << ' '
	 << setw(15) << "atoms" << endl;

    // loop over all trajectories
    int numFrames=0, readFrames=0;
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
       
        if (numFrames >= skip && (numFrames % stride) == 0) {
        (*pbc.*gathmethod)();

        double totSASA = 0;
        double totSASA_all = 0;

	for(unsigned int ir=0; ir < heavyatoms.size(); ++ir){

	  // determine the neighbours using a simple pairlist
	  // as a cutoff we use 2*rmax, because that is the furthest particles can be away and still 
	  // have their spheres overlap. It seems quite big, though. 

	  // Consider doing this not with the SimplePairlist, but with a loop over heavyatoms
	  // we could check that the distance is less than the sum of the radii (which should reduce 
	  // the number and we don't need to filter out the heavyatoms afterwards.
	  
	  
	  utils::AtomSpecifier pairlist(sys);
	  std::vector<int> neighbour;
	  for(unsigned int i=0; i < heavyatoms.size(); ++i){
	    if(i!=ir){
		
	      Vec dist=heavyatoms.pos(i) - heavyatoms.pos(ir);
	      if(dist.abs() < (radheavy[i] + radheavy[ir])){
		pairlist.addAtom(heavyatoms.mol(i), heavyatoms.atom(i));
		neighbour.push_back(i);
	      }
	    }
	  }
	  /*
	  utils::SimplePairlist pairlist(sys, *pbc, 2*rmax);
	  pairlist.setType("ATOMIC");
	  pairlist.setAtom(heavyatoms.mol(ir), heavyatoms.atom(ir));
	  pairlist.calc();

	  // we want to reduce the pairlist before we start looping over it (this is within 
	  // the loop over slices), to take only those into account that are also in the heavy atoms.
	  // we also keep track of the index of the atom in the heavyatoms 
	  utils::AtomSpecifier reduced_pairlist(sys);
	  std::vector<int> neighbour;
	  
	  for(int i = 0; i < pairlist.size(); ++i){
	      int n = heavyatoms.findAtom(pairlist.mol(i), pairlist.atom(i));
	      if(n>=0){
		reduced_pairlist.addAtom(pairlist.mol(i), pairlist.atom(i));
		neighbour.push_back(n);
	      }
	  }
	  */

          // set some variables

          double area = 0.0; // sums the area
          
          double rr = radheavy[ir];
          double rrsq = radheavysq[ir];

          if (pairlist.size() > 0) { 
	    // we have some neighbors so we start slicing up the atom

	    // determine number of slices to be done and the starting position
	    int nzp=(int) rint( (rr * 2) / zslice);
            double zgrid = heavyatoms.pos(ir)[2] - rr - zslice / 2;
            
            // section atom spheres perpendicular to the z axis
            // main inner loop
            for (int i = 0; i < nzp; ++i) {
              bool breakmain = false;
              double arcsum = 0; // sums the length of the arc
              zgrid += zslice;

              //     find the radius of the circle of intersection of 
              //     the ir sphere on the current z-plane
              double rsec2r = rrsq - (zgrid 
                     - heavyatoms.pos(ir)[2]) * (zgrid - heavyatoms.pos(ir)[2]);

	      // This should never be less than 0, but just to be sure
              double rsecr = 0.0;
	      if(rsec2r > 0) rsecr = sqrt(rsec2r);

              // vectors to store the start and end points of the arcs
              vector<double> arcf;
              vector<double> arci;

              // inner loop over neighbors
              for (unsigned int j = 0; j < pairlist.size(); ++j) {
                
		// calculate some distances
		Vec tmp = heavyatoms.pos(neighbour[j]);
		double dx = heavyatoms.pos(ir)[0] - tmp[0];
		double dy = heavyatoms.pos(ir)[1] - tmp[1];
		double dsq = dx * dx + dy * dy;
		double d = sqrt(dsq);
		
		//find radius of circle locus
		tmp = pairlist.pos(j);
		
		double rsec2n = radheavysq[neighbour[j]] - ((zgrid - tmp[2]) * (zgrid - tmp[2]));
		double rsecn=0.0;
		if(rsec2n > 0.0) rsecn = sqrt(rsec2n);
		  
		// find intersections of n.circles with ir circles in section
		// do the circles intersect, or is one circle completely inside the other?
		
		// this checks that 
		//  - at this z-coordinate we 'see the sphere': rsec2n > 0.0
		//  - the spheres can intersect: d < (rsecr + rsecn)
		if (rsec2n > 0.0 && d < (rsecr + rsecn)) {
		  
		  double diff_rsec = rsecr - rsecn;
		  
		  // this checks if 
		  //  - one circle is not completely within the other: d > abs(diff_rsec)
		  if (d > abs(diff_rsec)) {
		    
		    // find the points of intersection
		    
		    //     Initial and final arc endpoints are found for the ir circle intersected
		    //     by a neighboring circle contained in the same plane. The initial endpoint
		    //     of the enclosed arc is stored in arci, and the final arc in arcf
		    //     law of cosines
		    double trig_test = (dsq + rsec2r - rsec2n) / (2 * d * rsecr);

		    // If the spheres just touch (d=rsecr + rsecn), or are just within each 
		    // other (d=abs(diff_rsec)), this could be numerically more than 1 or less 
		    // than -1. Catch these. It should not really happen, though.
		    if (trig_test > 1.0) trig_test = 1.0;
		    if (trig_test < -1.0) trig_test = -1.0;

		    // alpha is the angle between a line containing a point of intersection and
		    // the reference circle center and the line containing both circle centers
		    double alpha = acos(trig_test);
		    
		    // beta is the angle between the line containing both circle centers and the x-axis
		    double beta = atan(dy/dx);

		    // atan gives a value between -PI/2 and +PI/2. We have to place it in the
		    // correct quadrant
		    if((dx < 0 && dy >= 0) || (dx <0 && dy < 0)) beta += PI;

		    // the arc that of the intersections goes from ti to tf 
		    double ti = beta - alpha;
		    double tf = beta + alpha;
		    
		    // keep them both between 0 and 2PI
		    if (ti < 0.0) ti += twoPI;
		    if (ti > twoPI) ti -= twoPI;
		    if (tf < 0.0)   tf += twoPI;
		    if (tf > twoPI) tf -= twoPI;
		    
		    arci.push_back(ti);
		    
		    if (tf < ti) {
		      //if the arc crosses zero, then it is broken into two segments.
		      //the first ends at twoPI and the second begins at zero
		      arcf.push_back(twoPI);
		      arci.push_back(0);
		    }
		    arcf.push_back(tf);
		    
		    
		  } else if(diff_rsec <=0) {
		    // - one circle is entirely in the other and ir is smaller than j:
		    //   this circle of ir does not contribute at all
		    j = pairlist.size();
		    breakmain = true;
		  }
		} // rsec2n > 0.0 && d < (rsecr + rsecn)
              } // j loop
	      
              //find the accessible surface area for the sphere ir on this section
	      // if breakmain == true: the sphere was completely within some other one
              if (!breakmain){
		
		// we have intersections, then we sum the contributions
		if(arci.size()){

		  //The arc endpoints are sorted on the value of the initial arc endpoint
		  vector<int> tag(arci.size(), 0);
		  heapsort(&arci[0], arci.size(), &tag[0]);
                  
		  arcsum = arci[0];
		  double t = arcf[tag[0]];
		  
		  for (unsigned int k = 1; k < arci.size(); ++k) {
		    if (t < arci[k]) arcsum += arci[k] - t;
		    double tt = arcf[tag[k]];
		    if (tt > t) t = tt;
		  }

		  arcsum += twoPI - t;

                } else { // no overlap with other circles, so it counts completely
                  arcsum = twoPI;
                }

		//calculate the accessible area
		
		// The area is equal to the accessible arc length x the radius of the 
		// circle x the section thickness.
		double parea = arcsum * rsecr * zslice;

                //Add the accessible area for this atom in this section to the area for this
                //atom for all the sections encountered thus far
                area += parea;
                
              } // breakmain = false
	    }// i loop
            
          } else { // we don't have neighbors, so calculate the exact area of the sphere

            area = 4*PI * rrsq;
          }

          // add it for averaging
          accs[ir] += area;
          if (atoms.findAtom(heavyatoms.mol(ir), heavyatoms.atom(ir)) >= 0) totSASA += area;
          totSASA_all += area;
          
        } // loop over atoms: ir
        
        cout.precision(5);
        cout << setw(10) << time << ' '
	     << setw(15) << totSASA << ' '
	     << setw(15) << totSASA_all << endl;
        readFrames++;
        } // skip/stride
        numFrames++;
      } // loop over frames in a single trajectory
      ic.close();
      
    } // loop over trajectories
    
    // calculate and print averages
    double totSASA = 0.0;
    double totSASA_all = 0.0;

    // loop over all heavy atoms
    for (unsigned int i = 0; i < heavyatoms.size(); ++i) {
      accs[i] /= readFrames;
      totSASA_all += accs[i];
      if (atoms.findAtom(heavyatoms.mol(i), heavyatoms.atom(i)) >= 0) totSASA += accs[i];
    }

    cout.precision(10);
    cout << "#\n# ave.          "
	 << setw(15) << totSASA << ' '
	 << setw(15) << totSASA_all << endl;
    if (args.count("verbose") >= 0) {

      cout.precision(5);
      cout << "#\n# average contribution per selected heavy atom\n";
      cout << "#\n# "
              << setw(6) << "atom"
              << setw(10) << "residue"
              << setw(5) << "name" << ' '
              << setw(15) << "SASA" << ' '
              << setw(15) << "\% selected" << ' '
              << setw(15) << "\% sasa atoms" << endl;

      double sumaccs = 0.0, sumper = 0.0, sumpera = 0.0;

      for (unsigned int i = 0; i < atoms.size(); ++i) {
        int index = heavyatoms.findAtom(atoms.mol(i), atoms.atom(i));
        if (index >= 0) {
          cout << "# "
                  << setw(6) << atoms.toString(i)
                  << setw(5) << atoms.resnum(i) + 1
                  << setw(5) << atoms.resname(i)
                  << setw(5) << atoms.name(i) << ' '
                  << setw(15) << accs[index] << ' '
                  << setw(15) << accs[index] / totSASA * 100.0 << ' '
                  << setw(15) << accs[index] / totSASA_all * 100.0 << ' '
                  << endl;
          sumaccs += accs[index];
          sumper += accs[index] / totSASA * 100.0;
          sumpera += accs[index] / totSASA_all * 100.0;
        }
      }
      cout << "#\n# total                 "
              << setw(15) << sumaccs << ' '
              << setw(15) << sumper << ' '
              << setw(15) << sumpera << ' ' << endl;
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void heapsort(double* values, int n, int* key) {

  //this implements a heapsort, which also returns the keys...
  //adapted to arrays that start indexing from 0...
  //     initialize index into the original ordering

  for (int i = 0; i < n; ++i) key[i] = i;
  //     perform the heapsort of the input list
  //		left = values.length/2;
  //		right = values.length-1;

  int k = n / 2 + 1;
  int index = n;
  double lists;
  int keys;

  do {

    if (k > 1) {
      k = k - 1;
      lists = values[k - 1];
      keys = key[k - 1];
    } else {
      lists = values[index - 1];
      keys = key[index - 1];
      values[index - 1] = values[0];
      key[index - 1] = key[0];
      index = index - 1;
    }
    if (index <= 1) {
      values[0] = lists;
      key[0] = keys;
      return;
    }

    int i = k;
    int j = k + k;
    do {

      if (j < index) {
        if (values[j - 1] < values[j]) ++j;
      }
      if (lists < values[j - 1]) {
        values[i - 1] = values[j - 1];
        key[i - 1] = key[j - 1];
        i = j;
        j = j + j;
      } else {
        j = index + 1;
      }
    } while (j <= index);

    values[i - 1] = lists;
    key[i - 1] = keys;

  } while (n > 1);


} //end void heapsort
