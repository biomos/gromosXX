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
 * @file jepot.cc
 * compute the 3J-value local elevation potential
 */

/**
 * @page programs Program Documentation
 *
 * @anchor jepot
 * @section jepot compute the 3J-value local elevation potential
 * @author @ref ja @ref mc
 * @date 01. 04. 09
 *
 * Program jepot computes the @f$^3J@f$-value local elevation (LE) potential from a LE
 * @f$^3J@f$-value restrained simulation. The LE potential can be calculated for all
 * values (@f$0 - 360^{\circ}@f$) of all restrained angles at the end of the simulation
 * only (\@fin) or for selected angles (\@angles) as a time-series throughout
 * the simulation (requires \@topo, \@pbc, \@postraj and \@restraj).
 * The \@timespec, \@timepts and \@restraj arguments control the time-series.
 * The time-series can be of the LE potential for all values of the selected angle
 * (ALL; default) or for only the current value of the selected angle (CURR) at each point
 * in time, giving only the current contribution of the LE potential to the overall
 * potential energy of the selected angle. With CURR, the \@jval file must contain the
 * @f$^3J@f$-value specifications for the selected angle only.
 *
 * \@K is the force constant given in the MD input file. Note that this is multiplied
 * by WJVR, the weight factor in the \@jval file, during the calculation of the LE potential
 * to give @f$k^{(j)}@f$.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@jval</td><td>&lt;jvalue restraint specifications&gt; </td></tr>
 * <tr><td> \@K</td><td>&lt;force constant&gt; </td></tr>
 * <tr><td> \@ngrid</td><td>&lt;number of grid points&gt; </td></tr>
 * <tr><td> [\@angles</td><td>&lt;angles over which to compute the LE potential: ALL (default) or CURR&gt;] </td></tr>
 * <tr><td> [\@fin</td><td>&lt;file containing final coordinates (if not time-series)&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time dt" (optional and only if time-series)&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints at which to compute the LE potential: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints at which to compute the LE potential (if time-series and \@timespec is EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@topo</td><td>&lt;molecular topology file (if CURR)&gt;] </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;boundary type (if CURR)&gt;] </td></tr>
 * <tr><td> [\@postraj</td><td>&lt;position trajectory file(s) (if CURR)&gt;] </td></tr>
 * <tr><td> [\@restraj</td><td>&lt;restraint trajectory file(s) (if time-series)&gt;] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 jepot
    @jval     ex.jval
    @K        0.01
    @ngrid    36
    @angles   CURR
    @time     0 0.01
    @timespec EVERY
    @timepts  100
    @topo     ex.top
    @pbc      r
    @postraj  ex.tr
    @restraj  ex.trs

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <sstream>
#include <math.h>
#include <cstdlib>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/utils/RestrTraj.h"
#include "../src/utils/groTime.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace utils;
using namespace gcore;
using namespace bound;
using namespace gio;

bool computeJepot(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done);

class karplus {
public:
  unsigned int m_i;
  unsigned int m_j;
  unsigned int m_k;
  unsigned int m_l;
  double weight;
  double j0;
  double delta;
  double A;
  double B;
  double C;

  karplus(int i, int j, int k, int l) {
    m_i = i;
    m_j = j;
    m_k = k;
    m_l = l;
  }

  karplus(const karplus &k) :
  m_i(k.m_i), m_j(k.m_j), m_k(k.m_k), m_l(k.m_l), weight(k.weight),
  j0(k.j0), delta(k.delta), A(k.A), B(k.B), C(k.C) {
  }
};

int main(int argc, char **argv) {

  Argument_List knowns;

  knowns << "jval" << "K" << "ngrid" << "angles" << "fin" << "time" << "timespec"
          << "timepts" << "topo" << "pbc" << "postraj" << "restraj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@jval        <jvalue restraint specifications>\n";
  usage += "\t@K           <force constant>\n";
  usage += "\t@ngrid       <number of grid points>\n";
  usage += "\t[@angles     <angle values over which to compute the LE potential: ALL (default) or CURR>]\n";
  usage += "\t[@fin         <file containing final coordinates (if not time-series)>]\n";
  usage += "\t[@time       <time dt (optional and only if time-series)>]\n";
  usage += "\t[@timespec   <timepoints at which to compute the LE potential: ALL (default), EVERY or SPEC (if time-series)>]\n";
  usage += "\t[@timepts    <timepoints at which to compute the LE potential (if time-series and timespec EVERY or SPEC)>]\n";
  usage += "\t[@topo       <molecular topology file (if CURR)>]\n";
  usage += "\t[@pbc        <boundary type (if CURR)>]\n";
  usage += "\t[@postraj    <position trajectory files (if CURR)>]\n";
  usage += "\t[@restraj    <restraint trajectory files (if time-series)>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // first check if we want the potential for all angle values at each point in time
    // or just for the current angle value
    string angles = "ALL";
    if (args.count("angles") > 0) {
      angles = args["angles"];
      if (angles != "ALL" && angles != "CURR")
        throw gromos::Exception("jepot",
              "angle specification " + angles + " unknown.\n");
    }

    // read in the j-value specifications
    Ginstream jf(args["jval"]);
    vector<string> buffer;
    jf.getblock(buffer);

    if (buffer[0] != "JVALRESSPEC")
      throw gromos::Exception("main", "jval file does not contain an JVALRESSPEC block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("jval", "J-value file " + jf.name() +
            " is corrupted. No END in " + buffer[0] +
            " block. Got\n"
            + buffer[buffer.size() - 1]);
    // at present can only handle one angle at a time if CURR
    if (angles == "CURR" && buffer.size() > 3)
      throw gromos::Exception("jepot", "To compute current-angle jepot, jval "
            "specification file should only contain one restraint\n");

    // store atom numbers, jvalues and karplus relation components
    // use the jvalue spec to create a property specifier for the torsion angle
    // because we need the angle spec to match the jval file anyway
    vector<karplus> kps;
    string spec = "t%1:";
    for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {
      // first get atoms as a string to define the property specifier
      if (angles == "CURR") {
        stringstream as(buffer[jj]);
        string ai, aj, ak, al;
        as >> ai >> aj >> ak >> al;
        spec = spec + ai + "," + aj + "," + ak + "," + al;
        if (as.fail())
          throw gromos::Exception("jepot", "Bad line in jval-file\n" + buffer[jj]);
      }

      // now get the data as integers, floats, etc
      istringstream is(buffer[jj]);
      int i, j, k, l;
      is >> i >> j >> k >> l;
      karplus kp(i, j, k, l);

      is >> kp.weight >> kp.j0 >> kp.delta >> kp.A >> kp.B >> kp.C;

      if (is.fail())
        throw gromos::Exception("jepot", "Bad line in jval-file\n" + buffer[jj]);
      kps.push_back(kp);
    }
    jf.close();

    // get K and n
    double K = args.getValue<double>("K");
    int ngrid = args.getValue<int>("ngrid");;
    // set bin width
    const double binwidth = 2.0 * M_PI / ngrid;

    // get time
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("jepot",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("jepot",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("jepot",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // check if we want a time-series (look for trajectory/ies)
    bool je_ts = false;
    if (args.count("restraj") > 0)
      je_ts = true;
    if (!je_ts) {
      if (args.count("fin") < 0) {
        throw gromos::Exception("jepot",
                "either specify @restraj or @fin:");
      } else if (angles == "CURR") {
        throw gromos::Exception("jepot",
                "CURR does not work with @fin:"
                " use @postraj and @restraj instead and use"
                " @timespec and @timepts to select the last frame");
      }
    }

    // if current angle contribution requested, define the properties
    // (torsion angles) based on the jval atoms
    // have to define things outside the if loop so they are accessible later in the program
    
    //PropertyContainer props(sys, pbc);
    vector<int> dihedral_types;
    vector<double> dihedral_angles;
    if (angles == "CURR") {

      //  read topology
      if (args.count("topo") > 0 && args.count("pbc") > 0 && args.count("postraj") > 0) {
        InTopology it(args["topo"]);
        System sys(it.system());
        GromosForceField gff(it.forceField());

        System refSys(it.system());

        // parse boundary conditions
        Boundary *pbc = BoundaryParser::boundary(sys, args);
        // parse gather method
        Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

        // define the properties
        PropertyContainer props(sys, pbc);
        {
          props.addSpecifier(spec);
          if (props[0]->type() != "Torsion")
            throw
            Arguments::Exception("Only dihedral (torsion) properties allowed");
        }

        // loop over all properties and store the dihedral types
        for (unsigned int i = 0; i < props.size(); i++) {
          int t = props[i]->getTopologyType(sys);
          if (t < 0)
            throw gromos::Exception("jepot", "Property " + props[i]->toTitle() +
                  " not found");
          dihedral_types.push_back(t);
        }

        // get the value of the dihedral angle at each point in time
        // set counters for angle timespecs
        int numTimepointsA = 0;
        // number of frames that have been written (actually number of angles calculated)
        unsigned int timesWrittenA = 0;
        // for SPEC: so that we stop trying when all requested timepoints are written
        bool doneA = false;

        // define input coordinate
        InG96 ic;
        // loop over all position trajectories
        for (Arguments::const_iterator
          iter = args.lower_bound("postraj"),
                to = args.upper_bound("postraj");
                iter != to; ++iter) {
          // open file    System sys;
    
          ic.open((iter->second).c_str());
          // loop over single position trajectory
          while (!ic.eof()) {
            ic >> sys;
            (*pbc.*gathmethod)();
            // update number for timespec
            numTimepointsA++;
            // check whether to use this block or move on
            if (computeJepot(numTimepointsA, timespec, timepts, timesWrittenA, doneA)) {
              // calculate the props
              props.calc();
              // store in the time series
              for (unsigned int i = 0; i < props.size(); i++) {
                dihedral_angles.push_back(props[i]->getValue().scalar());
              }
            }
            if (doneA)
              break;
          }
          ic.close();
        }
      } else {
        throw gromos::Exception("jepot", "To compute current-angle jepot, you must specify"
                " a topology file and the periodic boundary conditions\n");
      }
      // finally, initialise the output
      cout << setw(15) << "#time" << setw(18) << "phi" << setw(18) << "V_J_LE" << "\n";
    } // end if CURR


    // now compute the jepot for the time-series
    if (je_ts == true) {

      // set counters for jepot timespec
      int numTimepointsJ = 0;
      // number of frames that have been written
      unsigned int timesWrittenJ = 0;
      // for SPEC: so that we stop trying when all requested timepoints are written
      bool doneJ = false;

      // define input file type
      RestrTraj je;
      // loop over all restraint trajectories
      for (Arguments::const_iterator
        iter = args.lower_bound("restraj"),
              to = args.upper_bound("restraj");
              iter != to; ++iter) {

        // open this restraint trajectory file
        je.open((iter->second).c_str());
        // read in this restraint trajectory file
        while (!je.eof()) {
          je.read();
          // store jvalue eps data (for all restraints and time-points)
          JValueRestrData epsdata;
          je >> epsdata >> time;
          // update number for timespec
          numTimepointsJ++;
          // check whether to use this block or move on
          if (computeJepot(numTimepointsJ, timespec, timepts, timesWrittenJ, doneJ)) {

            // loop through the eps data
            for (unsigned int i = 0; i < epsdata.data().size(); ++i) {
              // check that it matches the kps data
              for (unsigned int j = 0; j < kps.size(); ++j) {
                // only take eps data that matches kps data (in case fewer restraints in kps)
                if ((epsdata.data()[i].i == kps[j].m_i) && (epsdata.data()[i].j == kps[j].m_j)
                        && (epsdata.data()[i].k == kps[j].m_k)
                        && (epsdata.data()[i].l == kps[j].m_l)) {

                  // if CURR only use one angle value
                  if (angles == "CURR") {
                    // convert dihedral angle to radians
                    double phi = (dihedral_angles[timesWrittenJ - 1] / 180.) * M_PI;
                    // compute contribution from each bin to this phi
                    double V = 0.0;
                    for (int bin = 0; bin < ngrid; ++bin) {
                      // phi0 = midpoint of this bin
                      const double phi0 = (bin + 0.5) * binwidth;
                      // correct for periodicity
                      double phi_bin = phi;
                      while (phi_bin < phi0 - M_PI)
                        phi_bin += 2 * M_PI;
                      while (phi_bin > phi0 + M_PI)
                        phi_bin -= 2 * M_PI;
                      // distance from corrected phi to midpoint of this bin
                      const double delta_phi = phi_bin - phi0;
                      // compute contribution from this phi to this bin
                      V += K * kps[j].weight * epsdata.data()[i].epsilon[bin] * exp(-delta_phi * delta_phi / (2 * binwidth * binwidth));
                    } // end loop over bins
                    // print contribution from all bins to the potential for this phi
                    cout << setw(18) << time << setw(18) << 180.0 * phi / M_PI << setw(18) << V << "\n";
                  } else {
                    // computing jepot for ALL angle values
                    // write time and atom numbers
                    cout << "#time:\t" << time << "\t\tatoms:\t" << epsdata.data()[i].i <<
                            "\t" << epsdata.data()[i].j << "\t" << epsdata.data()[i].k <<
                            "\t" << epsdata.data()[i].l << "\n\n";
                    // loop through possible phi values in 1deg increments (phi is in rads)
                    for (double phi = 0.0; phi < 2 * M_PI; phi += M_PI / 180.0) {
                      // compute contribution from each bin to this phi
                      double V = 0.0;
                      for (int bin = 0; bin < ngrid; ++bin) {
                        // phi0 = midpoint of this bin
                        const double phi0 = (bin + 0.5) * binwidth;
                        // correct for periodicity
                        double phi_bin = phi;
                        while (phi_bin < phi0 - M_PI)
                          phi_bin += 2 * M_PI;
                        while (phi_bin > phi0 + M_PI)
                          phi_bin -= 2 * M_PI;
                        // distance from corrected phi to midpoint of this bin
                        const double delta_phi = phi_bin - phi0;
                        // compute contribution from this phi to this bin
                        V += K * kps[j].weight * epsdata.data()[i].epsilon[bin] * exp(-delta_phi * delta_phi / (2 * binwidth * binwidth));
                      } // end loop over bins
                      // print contribution from all bins to the potential for this phi
                      cout << setw(18) << 180.0 * phi / M_PI << setw(18) << V << "\n";
                    } // end loop over phi
                    cout << "\n\n";
                  } // end else ALL
                } // end if atoms match loop
              } // end kps loop
            } // end epsdata loop
          } // end if for stride check
          if (doneJ)
            break;
        } // end time loop (eof loop)
      } // end loop over trajectories

      // read the epsilons in the case of only computing the final LE potential
      // we don't use the RestrTraj read function here because the format is different
    } else {
      // define input file
      Ginstream jpot(args["fin"]);
      // read in blocks until we find the JVALUERESEPS one
      jpot.getblock(buffer);
      while (buffer[0] != "JVALUERESEPS")
        jpot.getblock(buffer);

      if (buffer[buffer.size() - 1].find("END") != 0) {
        throw gromos::Exception("jepot", "Final coordinate file " + jpot.name() +
                " is corrupted. No END in " + buffer[0] +
                " block. Got\n"
                + buffer[buffer.size() - 1]);
      }

      // check if we have same number of jvalues in @jval file
      if (buffer.size() - 2 != kps.size()) {
        throw gromos::Exception("jepot", "J-value specification file has"
                " different number of restraints to final jepot file " + jpot.name());
      }

      // loop over jvalues
      for (unsigned int i = 1; i < buffer.size() - 1; i += 1) {
        // get epsilon for each bin
        istringstream is(buffer[i]);
        vector<double> epsilon(ngrid);
        for (int j = 0; j < ngrid; ++j)
          is >> epsilon[j];

        // write restraint number and atom numbers (read from jval file)
        cout << "\n\n# atoms: " << kps[i - 1].m_i << " "
                << kps[i - 1].m_j << " " << kps[i - 1].m_k << " " << kps[i - 1].m_l << " " << "\n";

        // loop through possible phi values in 1deg increments
        for (double phi = 0.0; phi < 2 * M_PI; phi += M_PI / 180.0) {
          // compute contribution from each bin to this phi
          double V = 0.0;
          for (int bin = 0; bin < ngrid; ++bin) {
            // phi0 = midpoint of this bin
            const double phi0 = (bin + 0.5) * binwidth;
            // correct for periodicity
            double phi_bin = phi;
            while (phi_bin < phi0 - M_PI)
              phi_bin += 2 * M_PI;
            while (phi_bin > phi0 + M_PI)
              phi_bin -= 2 * M_PI;
            // distance from corrected phi to midpoint of this bin
            const double delta_phi = phi_bin - phi0;
            // compute contribution from this phi to this bin
            V += K * kps[i - 1].weight * epsilon[bin] * exp(-delta_phi * delta_phi / (2 * binwidth * binwidth));
          } // end loop over bins
          // print contribution from all bins to the potential for this phi
          cout << setw(18) << 180.0 * phi / M_PI << setw(18) << V << "\n";
        } // end loop over phi
      } // end loop over j-value restraints
    } // end of if not timeseries

  } catch (const gromos::Exception &e) {
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool computeJepot(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done) {
  if (timespec == "ALL") {
    ++timesWritten;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesWritten;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesWritten;
        if (timesWritten == timepts.size())
          done = true;
        return true;
      } // compute jepot?
    } // times
  }
  return false;
}




