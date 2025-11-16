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
 * @file bilayer_gel.cc
 * Characterizes the geometry of a gel phase
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor bilayer_gel
 * @section bilayer_gel Characterizes the geometry of a gel phase
 * @author @ref bh
 * @date 28-10-09
 *
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@bilayeratoms</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  bilayer_gel
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @bilayeratoms     1-72:a
    @tail             12
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
#include <fstream>
#include <iomanip>
#include <math.h>

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
#include "../src/utils/VectorSpecifier.h"

#define PI 3.14159265

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
  knowns << "topo" << "pbc" << "time" << "bilayeratoms" << "tail"
          << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@bilayeratoms          <atomspecifier>\n";
  usage += "\t@tail          <integer: first atom of lipid tail>\n";
  usage += "\t@traj           <trajectory files>\n";
  
  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    // The two lines above read the topology and defines the system sys
    // sys is an instance of type System
    
    Boundary *pbc;
    pbc = BoundaryParser::boundary(sys, args);

    // get the @time argument
    utils::Time time(args);

    
    AtomSpecifier bilayer_atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("bilayeratoms");
      for(Arguments::const_iterator iter = args.lower_bound("bilayeratoms"); iter != to; iter++)
        bilayer_atoms.addSpecifier(iter->second);
    }
    
    int num_lipid = bilayer_atoms.mol(bilayer_atoms.size()-1) + 1;
    int half_lipid = num_lipid/2;
    vector<vector<Vec> > bigger(half_lipid);
    vector<vector<Vec> > smaller(half_lipid);
    int num_atperlip = (bilayer_atoms.size()+1)/num_lipid;

    cout << "Number of lipids: "<< num_lipid << endl;
    cout << "Number of atoms per lipid: "<< num_atperlip << endl;

    int tail = args.getValue<int>("tail", false, 0);
    if (tail < 0 || tail > num_atperlip)
      throw gromos::Exception(argv[0], "@tail: invalid range for first atom of tail.");

    utils::VectorSpecifier vs(sys,pbc);

    vector<double> times;
    // loop over all trajectories
    InG96 ic;
    int frames=0;
    Vec pointing(0.0,0.0,0.0);


    for(Arguments::const_iterator iter = args.lower_bound("traj"), to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while(!ic.eof()) {
        ic >> sys >> time;
        times.push_back(time.time());
        // gather system
        //(*pbc.*gathmethod)(); INBOX
        pbc->gather();
        Vec origin(sys.box().K()[0], sys.box().L()[1], sys.box().M()[2]);
        Vec shift = vs();
        origin /= 2;
        for(int i = 0; i < sys.numMolecules(); i++) {
          Vec cog(0.0, 0.0, 0.0);
          for(int j = 0; j < sys.mol(i).numAtoms(); j++) {
            sys.mol(i).pos(j) += shift;
            cog += sys.mol(i).pos(j);
          }
          cog /= sys.mol(i).numAtoms();
          cog = pbc->nearestImage(origin, cog, sys.box());
          for(int j = 0; j < sys.mol(i).numAtoms(); j++) {
            sys.mol(i).pos(j) = pbc->nearestImage(cog,
                    sys.mol(i).pos(j),
                    sys.box());
          }
        }
        for(int i = 0; i < sys.sol(0).numPos();
                i += sys.sol(0).topology().numAtoms()) {
          sys.sol(0).pos(i) += shift;

          sys.sol(0).pos(i) = pbc->nearestImage(origin,
                  sys.sol(0).pos(i),
                  sys.box());

          for(int j = 1; j < sys.sol(0).topology().numAtoms(); j++) {
            sys.sol(0).pos(i + j) += shift;

            sys.sol(0).pos(i + j) = pbc->nearestImage(sys.sol(0).pos(i),
                    sys.sol(0).pos(i + j),
                    sys.box());
          }
        }
        
        
        Vec cm = PositionUtils::com(sys,bilayer_atoms);
        int s = 0;
        int b = 0;
        for(int j=0;j<num_lipid;j++){
          Vec cml(0.0,0.0,0.0);
          double ml = 0;
          for(int a=0;a<num_atperlip;a++){
            int atom = (j*num_atperlip)+a;
            cml += bilayer_atoms.pos(atom) * bilayer_atoms.mass(atom);
            ml += bilayer_atoms.mass(atom);
          }
          int f_atom = (j*num_atperlip)+tail;
          int l_atom = (j*num_atperlip)+num_atperlip-1;
          Vec & btail = bilayer_atoms.pos(f_atom);
          Vec & etail = bilayer_atoms.pos(l_atom);
          pointing = btail - pbc->nearestImage(btail, etail, sys.box());
          pointing /=pointing.abs();
          cml /= ml;

          //cout << cm[2] << " " << cml[2] <<endl;

          //cout << pointing[0] << " " << pointing[1] << " "<<pointing[2]<<endl;

          if(cml[2] < cm[2]) {
            smaller[ s ].push_back(pointing);
            s++;
          } 
          else {
            bigger[ b ].push_back(pointing);
            b++;
          }
          //cout << smaller.size() << " " << bigger.size()<<endl;
        }
        frames++;

      } // while frames
      ic.close();
    } // for traj

    int counter_mono = 0;
    int counter_doub = 0;
    Vec av1(0.0,0.0,0.0);
    Vec av2(0.0,0.0,0.0);
    double phi_1=0;
    double phi_2=0;
    double theta_1=0;
    double theta_2=0;
    
    double sum_phi_1=0;
    double sum_phi_2=0;
    double sum_theta_1=0;
    double sum_theta_2=0;

    ofstream angles_ts;
    angles_ts.open("angles_ts.dat");
    angles_ts << "# Time series of the angles PHI1, PHI2, THETA1, THETA2\n"
            << "# " << endl
            << "# PHI1       PHI2       THETA1       THETA2\n";

    for(int it = 0; it < frames; it++) {// frames
      phi_1 = 0;
      phi_2 = 0;
      theta_1 = 0;
      theta_2 = 0;
      av1[0] = 0;
      av1[1] = 0;
      av1[2] = 0;
      av2[0] = 0;
      av2[1] = 0;
      av2[2] = 0;
      
      
      Vec unitz(0.0,0.0,1.0);
      counter_mono=0;
      counter_doub=0;
      for(int j = 0; j < half_lipid; j++) {// lipids
        av1+=bigger[j][it];
        av2+=smaller[j][it];
        counter_mono++;
      }
      av1/=av1.abs();
      av2/=av2.abs();
      

        phi_1 = atan2 (av1[1],av1[0]) * 180 / PI;
        phi_2 = atan2 (av2[1],av2[0]) * 180 / PI;

        if(phi_1<-80){
          phi_1 += 360;
        }
        if(phi_2<-80){
          phi_2 += 360;
        }

        //theta_1 = acos(av1[2]) * 180 / PI;

        theta_1 = acos(av1.dot(unitz)) * 180/PI;
        theta_2 = acos(av2.dot(-unitz)) * 180/PI;
        //theta2 += smaller[j][it].dot(-unitz);


       // prints components of pointing vectors
       //cout << av1[0] << "  " << av1[1] << "  " << av1[2] << endl;
       //cout << av2[0] << "  " << av2[1] << "  " << av2[2] << endl;
       


        

      sum_phi_1 += phi_1;
      sum_phi_2 += phi_2;
      sum_theta_1 += theta_1;
      sum_theta_2 += theta_2;

      angles_ts << times[it] << setw(14) << phi_1 << setw(14) << phi_2
              << setw(14) << theta_1 << setw(14) << theta_2 << endl;
    }
    //sum_phi /= (frames*counter_doub);
    sum_phi_1 /= (frames);
    sum_phi_2 /= (frames);
    sum_theta_1 /=(frames);
    sum_theta_2 /=(frames);
    
    
    cout << "AVERAGES:"<<endl;
    cout << "PHI1: "<< sum_phi_1 << endl;
    cout << "PHI2: "<< sum_phi_2 << endl;
    cout << "DeltaPHI: "<< sum_phi_1-sum_phi_2 << endl;
    cout << "THETA1: "<< sum_theta_1<< endl;
    cout << "THETA2: "<< sum_theta_2<< endl << endl;
    cout << "WARNING: Check the time series to verify jumping of the angles." << endl
            << " If there is jumping, the average values will be meaningless." << endl;


  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}

