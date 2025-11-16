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
 * @file ditrans.cc
 * monitors dihedral angle transitions
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ditrans
 * @section ditrans monitors dihedral angle transitions
 * @author @ref mc @ref co
 * @date 22. 11. 2004
 *
 * Dihedral angle transitions can be monitored during the course of a 
 * simulation using programs md and promd. Even though in many cases a 
 * molecular trajectory file will not contain every structure of the 
 * simulation, dihedral angle transitions can also be determined a posteriori
 * from a such a trajectory using program ditrans. This program can also write
 * the time series of dihedral angles without taking the inherent periodicity 
 * of a dihedral angle into account, but rather allow for dihedral angle values
 * below 0 degree or above 360 degree. 
 *
 * The program determines the position of maxima and minima in the dihedral 
 * angle potential energy function based on the phase shift and multiplicity 
 * given in the topology. Energy barriers arising from alternative terms, such
 * as nonbonded interactions, which may in theory shift the position of energy
 * minima and maxima of the dihedral angle are not taken into account. 
 *
 * Two different criteria can be used to count dihedral angle transitions, as 
 * described in the manual. A strict criterion only counts a transition once a
 * dihedral angle passes beyond the minimum of an adjacent energy well to 
 * prevent counting of short lived transitions of the maximum dividing the two
 * energy wells. Because of a possibly sparse sampling of data in a molecular
 * trajectory, this criterion may be too restrictive. As an alternative a
 * transition can also be counted as soon as a dihedral angle is seen to cross
 * the maximum separating two energy wells.
 *
 * The (standard) output can be restricted to the number of observed dihedral
 * angle transitions for every dihedral angle that was specified or can be
 * extended to information on every transition encountered.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref PropertySpecifier "properties"&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@strict</td><td>(use gromos96 transition criterion)] </td></tr>
 * <tr><td> [\@verbose</td><td>(print out every encountered transition)] </td></tr>
 * <tr><td> [\@tser</td><td>&lt;file name&gt; (extended time series)] </td></tr>
 * </table>
 *
 * <b>See also</b> @ref PropertySpecifier "property specifier"
 *
 * Example:
 * @verbatim
  ditrans
    @topo ex.top
    @pbc  r
    [@time 0 0.1]
    @prop t%1:1,3,5,6
    @traj ex.tr
 @endverbatim
 *
 * @bug Mar 22 2005: nearestImage calls in properties were missing
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Property.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Value.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

double nearest_minimum(double phi, double cosdelta, int multiplicity);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "verbose" << "prop" << "traj"
          << "strict" << "tser";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <boundary type>\n";
  usage += "\t@prop      <property specifier>\n";
  usage += "\t@traj      <trajectory files>\n";
  usage += "\t[@time      <T and dt>]\n";
  usage += "\t[@strict   (use gromos96 transition criterion)]\n";
  usage += "\t[@verbose  (print out every encountered transition)] \n";
  usage += "\t[@tser     <file name> (extended time series)]\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // check if verbose output required
    bool verb = false;
    if (args.count("verbose") >= 0)
      verb = true;

    // check if strict checking requested
    bool strict = false;
    if (args.count("strict") >= 0)
      strict = true;

    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);

    System sys(it.system());
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // get the @time argument
    utils::Time time(args);


    // read in a property

    // it's nice to actually store the properties read in somewhere
    // let's use a PropertyContainer, as it can handle some of the
    // work later on
    // the PropertyContainer should know about the system, so it can
    // check the input whether ie the atoms exist in the topology
    PropertyContainer props(sys, pbc);
    {
      Arguments::const_iterator iter = args.lower_bound("prop");
      Arguments::const_iterator to = args.upper_bound("prop");
      // we read in all properties specified by the user
      // all should be non-periodic dihedrals!!! (check that???)
      if (iter == to)
        throw Arguments::Exception("no property given");
      for (int i = 0; iter != to; iter++, ++i) {
        props.addSpecifier(iter->second);
        if (dynamic_cast<TorsionProperty*>(props[i]) == NULL)
          throw Arguments::Exception("Only dihedrals (non-periodic torsion; p) properties allowed");
      }
    }

    // do we want to write out the extended time series
    ofstream tser;
    bool do_tser = false;
    if (args.count("tser") > 0) {
      tser.open(args["tser"].c_str());
      do_tser = true;
      tser << "#" << setw(9) << "time";
      tser << "\t\t" << props.toTitle();
      tser << endl;
    }

    // loop over all properties and store the dihedral types
    vector<int> dihedral_types;
    vector<int> number_transitions;

    vector<double> old_min;
    vector<double> ts_old_val;
    vector<int> ts_offset;

    for (unsigned int i = 0; i < props.size(); i++) {
      int t = props[i]->getTopologyType(sys);
      if (t < 0)
        throw gromos::Exception("ditrans", "Property " + props[i]->toTitle() +
              " not found");
      dihedral_types.push_back(t);
      cout << "# property: " << props[i]->toTitle() << " type " << t+1 << endl;
      old_min.push_back(0.0);
      number_transitions.push_back(0);
      if (do_tser) {
        ts_offset.push_back(0);
        ts_old_val.push_back(0.0);
      }
    }

    // define input coordinate
    InG96 ic;
    bool lfirst = true;

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        (*pbc.*gathmethod)();

        // calculate the props
        // this is now the place, where a property-container is very handy
        // it knows that we want to loop over all properties and calculate
        // their 'value'. This works best if the property can be represented
        // by one number (which is often the case)
        props.calc();
        if (lfirst) {
          lfirst = false;
          for (unsigned int i = 0; i < props.size(); i++) {
            old_min[i] = nearest_minimum(props[i]->getValue().scalar(),
                    gff.dihedralType(dihedral_types[i]).pd(),
                    gff.dihedralType(dihedral_types[i]).np());
            if (do_tser) ts_old_val[i] = props[i]->getValue().scalar();
          }
        } else {
          if (!strict) {
            for (unsigned int i = 0; i < props.size(); i++) {
              double min = nearest_minimum(props[i]->getValue().scalar(),
                      gff.dihedralType(dihedral_types[i]).pd(),
                      gff.dihedralType(dihedral_types[i]).np());

              if (min != old_min[i]) {
                if (verb)
                  cout << "# at time " << time.time() << ": transition of "
                  << props[i]->toTitle() << " from " << old_min[i]
                        << " to " << min << endl;

                number_transitions[i]++;
                old_min[i] = min;

              }
            }
          } else {
            for (unsigned int i = 0; i < props.size(); i++) {

              double delta_phi = 360.0 / gff.dihedralType(dihedral_types[i]).np();
              double diff = utils::abs(old_min[i] - props[i]->getValue().scalar());
              if (diff > delta_phi && (360 - diff) > delta_phi) {

                //if (abs(old_min[i] - props[i]->getValue().scalar()) > delta_phi){
                double min =
                        nearest_minimum(props[i]->getValue().scalar(),
                        gff.dihedralType(dihedral_types[i]).pd(),
                        gff.dihedralType(dihedral_types[i]).np());

                if (verb)
                  cout << "# at time " << time.time() << ": transition of "
                  << props[i]->toTitle() << " from " << old_min[i]
                        << " to " << min << endl;

                number_transitions[i]++;
                old_min[i] = min;

              }
            }
          }
          // now possibly do the extended time series
          if (do_tser) {
            tser << time << "\t\t";
            for (unsigned int i = 0; i < props.size(); i++) {
              if (props[i]->getValue().scalar() < ts_old_val[i] - 180.0)
                ts_offset[i]++;
              if (props[i]->getValue().scalar() > ts_old_val[i] + 180.0)
                ts_offset[i]--;
              ts_old_val[i] = props[i]->getValue().scalar();
              tser << props[i]->getValue().scalar() + ts_offset[i]*360.0 << "\t\t";
            }
            tser << endl;
          }
        }
      }

      ic.close();
    }
    tser.close();
    for (unsigned int i = 0; i < props.size(); i++) {
      cout << props[i]->toTitle() << "\t" << number_transitions[i] << endl;
    }

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

double nearest_minimum(double phi, double cosdelta, int multiplicity) {
  double a_minimum = 180.0 * (3.0 - cosdelta) / (2.0 * multiplicity);
  double delta_phi = 360.0 / multiplicity;
  //  double minimum_within_pi=a_minimum - int(rint((a_minimum - phi)/delta_phi))*2*M_PI;
  double nearest_min = a_minimum
          - int(rint((a_minimum - phi) / delta_phi)) * delta_phi;
  // not so nice to get it down to zero again fi cosdelta = -1 and nearest_min
  // is almost 360.0
  if (nearest_min >= 360.0 - 0.001) nearest_min -= 360.0;

  return nearest_min;
}

