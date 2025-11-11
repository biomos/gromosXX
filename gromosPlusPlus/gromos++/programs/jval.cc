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
 * @file jval.cc
 * compute 3J-values and statistics
 */

/**
 * @page programs Program Documentation
 *
 * @anchor jval
 * @section jval compute 3J-values and statistics
 * @author @ref co @ref mc @ref ja
 * @date 22. 11. 2004
 *
 * Program jval calculates @f$^3J@f$-values from a single conformation or from a
 * trajectory. It can write out the values of all @f$^3J@f$-couplings
 * specified in the \@jval file or the total rmsd over all couplings
 * from the reference values at each point in time. The final part of the output
 * is always a summary of the @f$^3J@f$-value specification parameters, the
 * averages over the entire trajectory and other statistics.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type [&lt;gather method&gt;]&gt; </td></tr>
 * <tr><td> \@jval</td><td>&lt;@f$^3J@f$-value specification file&gt; </td></tr>
 * <tr><td> [\@timeseries</td><td>&lt;write time-series of @f$^3J@f$-values]&gt; </td></tr>
 * <tr><td> [\@rmsd</td><td>&lt;write the rmsd over all @f$^3J@f$-values as a time-series&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time dt" (optional and only if time-series)&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints at which to compute the @f$^3J@f$-values: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints at which to compute the @f$^3J@f$-values (if time-series and \@timespec is EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;position trajectory file(s)&gt; </td></tr>
 * </table>
 * 
 * Example:
 * @verbatim
   jval
     @topo ex.top
     @pbc  r
     @jval ex.jval
     @traj ex.tr
   @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Stat.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

// function to skip time-points
bool computeJval(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done);

// structure to store the information needed to define the Karplus relation
class karplus{
 public:
  int m_mol;
  int m_i;
  int m_j;
  int m_k;
  int m_l;
  double j0;
  double delta;
  double A;
  double B;
  double C;
  karplus(const System &sys, int i, int j, int k, int l)
    {
      //determine molecule
      int m=0, offset=0;
      
      for(; i>=sys.mol(m).numAtoms()+offset; m++)
	offset += sys.mol(m).numAtoms();
      m_mol=m;
      
      m_i=i-offset;
      m_j=j-offset;
      m_k=k-offset;
      m_l=l-offset;
    }
  karplus(const karplus &k):
    m_mol(k.m_mol), m_i(k.m_i), m_j(k.m_j), m_k(k.m_k), m_l(k.m_l), j0(k.j0),
    delta(k.delta), A(k.A), B(k.B), C(k.C) {}
};


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "jval" << "timeseries" << "rmsd" << "time" << "timespec" << "timepts" << "avg" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n";
  usage += "\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type [<gather method>]>\n";
  usage += "\t@jval         <j-value specification file>\n";
  usage += "\t[@timeseries  <write time-series of j-values>\n";
  usage += "\t[@rmsd        <write the rmsd over all j-values as a time-series>\n";
  usage += "\t[@time        <time and dt> (optional and only if time-series)]\n";
  usage += "\t[@timespec    <timepoints at which to compute the j-values: ALL (default), EVERY or SPEC (if time-series)]\n";
  usage += "\t[@timepts     <timepoints at which to compute the j-values (if time-series and timespec EVERY or SPEC)]\n";
  usage += "\t@traj         <position trajectory file(s)>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);
    
    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);

    System sys(it.system());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

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
    // to store karplus relation information
    vector<karplus> kps;

    for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {
      istringstream is(buffer[jj]);
      int i, j, k, l;
      double fdum;
      is >> i >> j >> k >> l >> fdum;
      karplus kp(sys, i - 1, j - 1, k - 1, l - 1);
      // cout << i << " " << j << " " << k << " " << l << " " << fdum;

      is >> kp.j0 >> kp.delta >> kp.A >> kp.B >> kp.C;
      //cout << kp.m_i << " " << kp.m_j << " " << kp.m_k << " " << kp.m_l << 
      //	" " << kp.j0 << " " << kp.delta << " " << kp.A << " " << kp.B 
      //   << " " << kp.C << " " << endl;

      if (is.fail())
        throw gromos::Exception("jval", "Bad line in jval-file\n" + buffer[jj]);
      kps.push_back(kp);
    }
    jf.close();

    // now we have to create the properties
    PropertyContainer props(sys, pbc);
    for (unsigned int i = 0; i < kps.size(); i++) {
      ostringstream os;
      os << "t%" << kps[i].m_mol + 1 << ":"
              << kps[i].m_i + 1 << "," << kps[i].m_j + 1 << ","
              << kps[i].m_k + 1 << "," << kps[i].m_l + 1;
      props.addSpecifier(os.str());
    }

    //   get simulation time if given
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("jval",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("jval",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("jval",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // check if we want a time-series of the j-values
    bool jv_ts = false;
    if (args.count("timeseries") >= 0)
      jv_ts = true;

    // check if we want to do a time-series of the rmsd
    bool rmsd_ts = false;
    if (args.count("rmsd") >= 0) {
      if (jv_ts) {
        throw gromos::Exception("jval",
                "cannot write a time-series of the j-values and the rmsd:"
                " use either @rmsd or @timeseries");
      } else {
        rmsd_ts = true;
        // write a heading
        cout << setw(8) << "#time" << setw(15) << "rmsd" << endl;
      }
    }


    /* old: for printing rmsd time-series
  bool avg_ts = false;
  if (args.count("avg") >= 0)
    avg_ts = true;
     */

    // set counters for timespecs
    int numTimepoints = 0;
    // number of frames that have been written
    unsigned int timesWritten = 0;
    // for SPEC: so that we stop trying when all requested timepoints are written
    bool done = false;

    // define input coordinate
    InG96 ic;

    // define enough statistic classes
    vector<gmath::Stat<double> > stat;
    stat.resize(2*kps.size());

    // loop over all trajectories
    for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
        /* old time stuff
        if (time.doSeries())*/
        ic >> time;
         
        (*pbc.*gathmethod)();

        // update number for timespec
        numTimepoints++;
        // check whether to use this block or move on
        if (computeJval(numTimepoints, timespec, timepts, timesWritten, done)) {

          // calculate the props
          props.calc();

          // this is for if we write out for each jval vs time
          /* old
                    if (time.doSeries() && !avg_ts) {
                      cout << "#\n#\n# TIME\t" << time << endl;
                      cout << "#" << setw(4) << "num"
                              << setw(12) << "j" << endl;
                    }
           */

          // title for printing jvalues at this time-point
          if (jv_ts && !rmsd_ts) {
            cout << "#\n#\n# TIME\t" << time << endl;
            cout << "#" << setw(4) << "num" << setw(16) << "j-value" << endl;
          }

          // calculate the j-values and store

          double rmsd = 0.0;

          for (unsigned int i = 0; i < kps.size(); i++) {
            double cosphi = cos((props[i]->getValue().scalar() + kps[i].delta) * M_PI / 180.0);
            double J = kps[i].A * cosphi * cosphi +
                    kps[i].B * cosphi +
                    kps[i].C;
            stat[i].addval(J);
            stat[kps.size() + i].addval(props[i]->getValue().scalar());

            //if (time.doSeries() && !avg_ts) { // old
            
            // write out each number and j-value
            if (jv_ts && !rmsd_ts) {
              cout << setw(5) << i + 1 << setw(16) << J << endl;
            }
            rmsd += (J - kps[i].j0) * (J - kps[i].j0);
          }

          rmsd /= kps.size();

          // write out time-series of rmsd
          //if (time.doSeries() && avg_ts) { // old
          if (rmsd_ts) {
            cout << setw(8) << time << setw(15) << sqrt(rmsd) << endl;
          }
        } // end if computeJval
        if (done)
          break;
      }
      ic.close();
    } // end loop over trajectories

    // now write out the averages (whether time-series or not)
    cout << "#\n#\n# JVALUE AVERAGES\n#\n"
            << "# (print with a2ps --landscape --columns=1 --font-size=9)\n#\n";

    // print title
    cout << "#"
            << setw(4) << "num"
            << setw(4) << "mol"
            << setw(9) << "residue"
            << setw(15) << "atom names"
            << setw(21) << "atom numbers"
            << setw(11) << "A"
            << setw(5) << "B"
            << setw(5) << "C"
            << setw(9) << "delta"
            << setw(7) << "J0"
            << setw(10) << "phi ave"
            << setw(9) << "phi rmsd"
            << setw(9) << "J ave"
            << setw(8) << "J rmsd"
            << setw(10) << "|Jave-J0|"
            << endl;

    // now print out and calculate overall performance
    double sum = 0;
    double abssum = 0;
    double ssum = 0;

    for (unsigned int i = 0; i < kps.size(); i++) {
      int m = kps[i].m_mol;
      cout << setw(5) << i + 1
              << setw(4) << m + 1
              << setw(4) << sys.mol(m).topology().resNum(kps[i].m_j) + 1
              << setw(5) << sys.mol(m).topology().resName(sys.mol(m).topology().resNum(kps[i].m_j))
              << setw(5) << sys.mol(m).topology().atom(kps[i].m_i).name()
              << setw(5) << sys.mol(m).topology().atom(kps[i].m_j).name()
              << setw(5) << sys.mol(m).topology().atom(kps[i].m_k).name()
              << setw(5) << sys.mol(m).topology().atom(kps[i].m_l).name()
              << setw(7) << kps[i].m_i + 1
              << setw(5) << kps[i].m_j + 1
              << setw(5) << kps[i].m_k + 1
              << setw(5) << kps[i].m_l + 1;

      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(1);

      cout << setw(7) << kps[i].A
              << setw(5) << kps[i].B
              << setw(5) << kps[i].C;
      cout.precision(0);
      cout << setw(7) << kps[i].delta;
      cout.precision(2);
      cout << setw(7) << kps[i].j0;
      cout.precision(1);
      cout << setw(10) << stat[kps.size() + i].ave();
      cout.precision(2);
      cout << setw(9) << stat[kps.size() + i].rmsd()
              << setw(9) << stat[i].ave();
      cout.precision(3);
      cout << setw(8) << stat[i].rmsd();
      cout.precision(1);
      cout << setw(10) << fabs(stat[i].ave() - kps[i].j0)
              << endl;
      sum += stat[i].ave() - kps[i].j0;
      abssum += fabs(stat[i].ave() - kps[i].j0);
      ssum += (kps[i].j0 - stat[i].ave())*(kps[i].j0 - stat[i].ave());
    }
    cout << "\n#"
            << setw(30) << "average deviation " << sum / kps.size() << endl;
    cout << "#"
            << setw(30) << "average absolute deviation "
            << abssum / kps.size() << endl;
    cout << "#"
            << setw(30) << "root-mean-square deviation "
            << sqrt(ssum / kps.size()) << endl;
    cout << "#"
            << setw(30) << "rmsd over deviations "
            << sqrt((ssum - sum * sum / kps.size()) / kps.size()) << endl;

  }  catch (const gromos::Exception &e) {
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool computeJval(int i, std::string const & timespec, vector<int> const & timepts,
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
      } // compute jval?
    } // times
  }
  return false;
}



