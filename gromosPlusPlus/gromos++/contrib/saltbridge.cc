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
 * @file saltbridge.cc
 * Monitors the occurrence of saltbridges
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor saltbridge
 * @section saltbridge Monitors the occurrence of saltbridges
 * @author @ref bh
 * @date 8-7-09
 *
 * Monitors the occurrence of saltbridges over a molecular
 * trajectory file through geometric criteria (cutoff distance between positive
 * and negative residues).
 *
 * For each type of saltbridge, a cutoff correction is stablished based on 
 * simple combination rules. For now, the only acceptable types are ARG-GLU,
 * ARG-ASP, LYSH-GLU, LYSH-ASP, HISH-GLU and HISH-ASP. These types are the most
 * common in protein-protein interactions. Correction factors were estimated
 * from the analysis of PDB structures. These values are yet hard coded.
 *
 * The user can specify two groups of atoms (donors and acceptors) between which
 * the saltbridge interactions are to be monitored. The program can be extended
 * to non-protein systems. However, this requires the implementation of a 
 * library file in order to define residue pairs and correction factors.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@donor</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@acceptor</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> [\@cutoff</td><td>&lt;distance [nm]; default: 0.5;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  saltbridge
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @donor            1:a
    @acceptor         2:a
    @cutoff            0.5 
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
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "donor" << "acceptor" << "cutoff"
          << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@donor          <atomspecifier>\n";
  usage += "\t@acceptor       <atomspecifier>\n";
  usage += "\t[@cutoff        <distance [nm]; default: 0.5>]\n";
  usage += "\t@traj           <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    // The two lines above read the topology and defines the system sys
    // sys is an instance of type System

    System refSys(it.system());

    // The two lines below define pbc as a pointer of type Boundary
    // and the gethering method
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // get the @time argument
    utils::Time time(args);

    // The donor_atoms is an instance of type AtomSpecifier.
    // It is for example one protein considered "donor"
    // Actually it doesn't matter who is donor and who is acceptor in this code
    AtomSpecifier donor_atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("donor");
      for (Arguments::const_iterator iter = args.lower_bound("donor"); iter != to; iter++)
        donor_atoms.addSpecifier(iter->second);
    }

    // Here is the same but now for acceptor
    // Donor and acceptor can be the same
    AtomSpecifier acceptor_atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("acceptor");
      for (Arguments::const_iterator iter = args.lower_bound("acceptor"); iter != to; iter++)
        acceptor_atoms.addSpecifier(iter->second);
    }

    // Now, for each group (donor or acceptor), positive and negative residues
    // are defined and collected. For example, array donor_pos will contain the
    // positive residues of donor.
    AtomSpecifier donor_pos(sys), donor_neg(sys);
    for (int i = 0; i < donor_atoms.size(); i++) {
      if (donor_atoms.resname(i) == "LYSH" && donor_atoms.name(i) == "NZ")
        donor_pos.addAtom(donor_atoms.mol(i), donor_atoms.atom(i));
      if (donor_atoms.resname(i) == "ARG" && donor_atoms.name(i) == "CZ")
        donor_pos.addAtom(donor_atoms.mol(i), donor_atoms.atom(i));
      if (donor_atoms.resname(i) == "HISH" && donor_atoms.name(i) == "CE1")
        donor_pos.addAtom(donor_atoms.mol(i), donor_atoms.atom(i));
      if (donor_atoms.resname(i) == "ASP" && donor_atoms.name(i) == "CG")
        donor_neg.addAtom(donor_atoms.mol(i), donor_atoms.atom(i));
      if (donor_atoms.resname(i) == "GLU" && donor_atoms.name(i) == "CD")
        donor_neg.addAtom(donor_atoms.mol(i), donor_atoms.atom(i));
    }


    AtomSpecifier acceptor_pos(sys), acceptor_neg(sys);
    for (int i = 0; i < acceptor_atoms.size(); i++) {
      if (acceptor_atoms.resname(i) == "LYSH" && acceptor_atoms.name(i) == "NZ")
        acceptor_pos.addAtom(acceptor_atoms.mol(i), acceptor_atoms.atom(i));
      if (acceptor_atoms.resname(i) == "ARG" && acceptor_atoms.name(i) == "CZ")
        acceptor_pos.addAtom(acceptor_atoms.mol(i), acceptor_atoms.atom(i));
      if (acceptor_atoms.resname(i) == "HISH" && acceptor_atoms.name(i) == "CE1")
        acceptor_pos.addAtom(acceptor_atoms.mol(i), acceptor_atoms.atom(i));
      if (acceptor_atoms.resname(i) == "ASP" && acceptor_atoms.name(i) == "CG")
        acceptor_neg.addAtom(acceptor_atoms.mol(i), acceptor_atoms.atom(i));
      if (acceptor_atoms.resname(i) == "GLU" && acceptor_atoms.name(i) == "CD")
        acceptor_neg.addAtom(acceptor_atoms.mol(i), acceptor_atoms.atom(i));
    }

    // get the cutoff parameter
    double cutoff_dis = args.getValue<double>("cutoff", false, 0.5);

    // Defining the corrections for salt-bridge type
    // The negative residues have the same functional group (carboxylate)
    // But the positive residues have different functional groups
    // and the distances are defined in a different way (reference atoms)
    // Thus, below are the distance corrections for cutoff
    // In the future, this could be included in a library file
    double lysh_corr = 0.115;
    double arg_corr = 0.05;
    double hish_corr = 0.02;
    double corr = 0; // this will take one of the above values

    // Below are the initialization of the distance matrices    
    vector<vector<double> > occupancy_matrix1(donor_pos.size(),
            vector<double>(acceptor_neg.size(), 0.0));

    for (int i = 0; i < donor_pos.size(); i++) {
      for (int j = 0; j < acceptor_neg.size(); j++) {
        occupancy_matrix1[i][j] = 0.0;
      }
    }

    vector<vector<double> > occupancy_matrix2(donor_neg.size(),
            vector<double>(acceptor_pos.size(), 0.0));

    for (int i = 0; i < donor_neg.size(); i++) {
      for (int j = 0; j < acceptor_pos.size(); j++) {
        occupancy_matrix2[i][j] = 0.0;
      }
    }

    int counterframe = 0;
    // loop over all trajectories
    InG96 ic;

    // create the ouput file for the time series
    ofstream timeseries("saltbridges_ts.dat");
    if (!timeseries.good())
      throw gromos::Exception("saltbridge", "Cannot write time series.");
    timeseries << "# Columns are: " << endl
            << "#   1: time" << endl;
    {
      int column_number = 2;
      for (int i = 0; i < donor_pos.size(); i++) {
        for (int j = 0; j < acceptor_neg.size(); j++, column_number++) {
          timeseries << "# " << setw(3) << column_number << ": "
                  << setw(5) << donor_pos.mol(i) + 1 << setw(5) << donor_pos.resnum(i) + 1
                  << setw(5) << donor_pos.resname(i) << " - "
                  << setw(5) << acceptor_neg.mol(j) + 1 << setw(5) << acceptor_neg.resnum(j) + 1
                  << setw(5) << acceptor_neg.resname(j) << endl;
        }
      }
      for (int i = 0; i < donor_neg.size(); i++) {
        for (int j = 0; j < acceptor_pos.size(); j++, column_number++) {
          timeseries << "# " << setw(3) << column_number << ": "
                  << setw(5) << donor_neg.mol(i) + 1 << setw(5) << donor_neg.resnum(i) + 1
                  << setw(5) << donor_neg.resname(i) << " - "
                  << setw(5) << acceptor_pos.mol(j) + 1 << setw(5) << acceptor_pos.resnum(j) + 1
                  << setw(5) << acceptor_pos.resname(j) << endl;
        }
      }
    }


    for (Arguments::const_iterator iter = args.lower_bound("traj"), to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        // gather system
        (*pbc.*gathmethod)();

        // Calculating and printing distances
        timeseries << time;
        for (int i = 0; i < donor_pos.size(); i++) {
          for (int j = 0; j < acceptor_neg.size(); j++) {
            Vec & pos_i = donor_pos.pos(i);
            Vec & pos_j = acceptor_neg.pos(j);
            Vec dist = pos_i - pbc->nearestImage(pos_i, pos_j, sys.box());
            double d = dist.abs();

            if (d < cutoff_dis) {
              if (donor_pos.resname(i) == "LYSH")
                corr = lysh_corr;
              if (donor_pos.resname(i) == "ARG")
                corr = arg_corr;
              if (donor_pos.resname(i) == "HISH")
                corr = hish_corr;
              double effect_dist = d + corr;
              if (effect_dist < cutoff_dis)
                occupancy_matrix1[i][j] += 1;
            }
            timeseries << setw(15) << d;
          }
        }

        for (int i = 0; i < donor_neg.size(); i++) {
          for (int j = 0; j < acceptor_pos.size(); j++) {
            Vec & pos_i = donor_neg.pos(i);
            Vec & pos_j = acceptor_pos.pos(j);
            Vec dist = pos_i - pbc->nearestImage(pos_i, pos_j, sys.box());
            double d = dist.abs();
            if (d < cutoff_dis) {
              if (acceptor_pos.resname(j) == "LYSH")
                corr = lysh_corr;
              if (acceptor_pos.resname(j) == "ARG")
                corr = arg_corr;
              if (acceptor_pos.resname(j) == "HISH")
                corr = hish_corr;
              double effect_dist = d + corr;
              if (effect_dist < cutoff_dis)
                occupancy_matrix2[i][j] += 1;
            }
            timeseries << setw(15) << d;
          }
        }
        timeseries << endl;


        counterframe += 1;
      } // while frames
      ic.close();
    } // for traj

    double min_occ = 0.1;
    cout << "# List of saltbridges with occupancy bigger than: "
            << min_occ << endl;
    cout << "# Matrix 1:" << endl;


    for (int i = -1; i < donor_pos.size(); i++) {
      for (int j = -1; j < acceptor_neg.size(); j++) {

        if (i != -1 && j != -1) {
          double occ_fraction = occupancy_matrix1[i][j] / counterframe;

          if (occ_fraction >= min_occ) {
            cout << "Molec:" << setw(3) << donor_pos.mol(i) + 1
                    << " " << setw(4) << donor_pos.resname(i)
                    << setw(3) << donor_pos.resnum(i) + 1
                    << " X "
                    << "Molec:" << setw(3) << acceptor_neg.mol(j) + 1
                    << " " << setw(4) << acceptor_neg.resname(j)
                    << setw(3) << acceptor_neg.resnum(j) + 1
                    << " :" << setw(15) << "Occupancy: " << occ_fraction << endl;
          }
        }
      }
    }

    cout << "# Matrix 2:" << endl;

    for (int i = -1; i < donor_neg.size(); i++) {
      for (int j = -1; j < acceptor_pos.size(); j++) {

        if (i != -1 && j != -1) {
          double occ_fraction = occupancy_matrix2[i][j] / counterframe;

          if (occ_fraction >= min_occ) {
            cout << "Molec:" << setw(3) << donor_neg.mol(i) + 1
                    << " " << setw(4) << donor_neg.resname(i)
                    << setw(3) << donor_neg.resnum(i) + 1
                    << " X "
                    << "Molec:" << setw(3) << acceptor_pos.mol(j) + 1
                    << " " << setw(4) << acceptor_pos.resname(j)
                    << setw(3) << acceptor_pos.resnum(j) + 1
                    << " :" << setw(15) << "Occupancy: " << occ_fraction << endl;
          }
        }
      }
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}

