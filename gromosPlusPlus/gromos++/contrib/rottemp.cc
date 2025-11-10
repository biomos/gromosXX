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
 * @file rottemp.cc
 * calculates temperature distributions
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor rottemp
 * @section rottemp calculates temperature distributions
 * @author unkown
 * @date unkown
 *
 * Program rottemp calculates internal, rotational and translational temperature
 * distributions.
 *
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions, gathering method&gt; </td></tr>
 * <tr><td> \@timestep</td><td>&lt;simulation timestep&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;position or combined trajectory files&gt;</td></tr>
 * <tr><td>[\@trav</td><td>&lt;velocity trajectory files&gt;]</td></tr>
 * <tr><td>[\@constr</td><td>&lt;switch constraints on/off&gt;]</td></tr>
 * </table>
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
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
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Distribution.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "timestep" << "traj" << "trav"
          << "constr";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <topology>\n";
  usage += "\t@pbc         <boundary type>\n";
  usage += "\t@timestep    <timestep of the simulation>\n";
  usage += "\t@traj        <position or combined trajectory files>\n";
  usage += "\t[@trav       <velocity trajectory files>]\n";
  usage += "\t[@constr     <switch constraints on / off>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);

    System sys(it.system());

    double timestep = args.getValue<double>("timestep", true);

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // position or combined trajectory
    InG96 itrj;

    // velocity trajectory
    InG96 itrv;
    bool read_trv = false;
    Arguments::const_iterator trv_iter, trv_to;
    if (args.count("trav") > 0) {
      read_trv = true;
      trv_iter = args.lower_bound("trav");
      trv_to = args.upper_bound("trav");
    }

    bool constr = false;
    if (args.count("constr") >= 0)
      constr = true;

    Stat<double> trans_temp, rot_temp, int_temp;

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      itrj.open((iter->second).c_str());

      if (read_trv) {
        assert(trv_iter != trv_to);
        itrv.open((trv_iter->second).c_str());
        ++trv_iter;
      }

      // loop over single trajectory
      while (!itrj.eof()) {
        itrj >> sys;
        // time seems not to be used
        //itrj >> sys >> time;
        // gather
        (*pbc.*gathmethod)();

        if (read_trv)
          itrv >> sys;
        // time seems not to be used
        //itrv >> sys >> time;


        for (int i = 0; i < sys.numMolecules(); ++i) {

          Vec com_v(0.0, 0.0, 0.0),
                  com_r(0.0, 0.0, 0.0);
          double com_mass = 0.0;
          double tot_ekin = 0.0;

          for (int a = 0; a < sys.mol(i).numAtoms(); ++a) {
            const double mass = sys.mol(i).topology().atom(a).mass();
            com_mass += mass;
            com_v += mass * sys.mol(i).vel(a);
            tot_ekin += 0.5 * mass * sys.mol(i).vel(a).dot(sys.mol(i).vel(a));
            com_r += mass * sys.mol(i).pos(a) -
                    0.5 * timestep * mass * sys.mol(i).vel(a);
          }

          com_v /= com_mass;
          com_r /= com_mass;

          // check...
          cout.precision(4);
          cout.setf(ios::fixed, ios::floatfield);

          double const ekin_trans = 0.5 * com_mass * com_v.dot(com_v);

          Vec com_L(0.0, 0.0, 0.0);
          Matrix com_I(3, 3, 0.0);
          Vec r;

          for (int a = 0; a < sys.mol(i).numAtoms(); ++a) {
            const double mass = sys.mol(i).topology().atom(a).mass();
            r = sys.mol(i).pos(a) -
                    0.5 * timestep * sys.mol(i).vel(a) - com_r;

            com_L += mass * sys.mol(i).pos(a).cross(sys.mol(i).vel(a));

            com_I(0, 0) += mass * (r[1] * r[1] + r[2] * r[2]);
            com_I(1, 1) += mass * (r[0] * r[0] + r[2] * r[2]);
            com_I(2, 2) += mass * (r[0] * r[0] + r[1] * r[1]);
            com_I(1, 0) += mass * (-r[0] * r[1]);
            com_I(0, 1) += mass * (-r[0] * r[1]);
            com_I(2, 0) += mass * (-r[0] * r[2]);
            com_I(0, 2) += mass * (-r[0] * r[2]);
            com_I(2, 1) += mass * (-r[1] * r[2]);
            com_I(1, 2) += mass * (-r[1] * r[2]);
          }

          com_L -= com_mass * com_r.cross(com_v);

          // invert the inertia tensor
          Matrix com_II(3, 3, 0);
          const double denom = -com_I(2, 0) * com_I(2, 0) * com_I(1, 1)
                  + 2 * com_I(0, 1) * com_I(0, 2) * com_I(1, 2)
                  - com_I(0, 0) * com_I(1, 2) * com_I(1, 2)
                  - com_I(0, 1) * com_I(0, 1) * com_I(2, 2)
                  + com_I(0, 0) * com_I(1, 1) * com_I(2, 2);

          com_II(0, 0) = (-com_I(1, 2) * com_I(1, 2) + com_I(1, 1) * com_I(2, 2));
          com_II(1, 0) = com_II(0, 1) = (com_I(0, 2) * com_I(1, 2)
                  - com_I(0, 1) * com_I(2, 2));
          com_II(0, 2) = com_II(2, 0) = (-com_I(0, 2) * com_I(1, 1)
                  + com_I(0, 1) * com_I(1, 2));

          com_II(1, 1) = (-com_I(0, 2) * com_I(0, 2) + com_I(0, 0) * com_I(2, 2));
          com_II(1, 2) = com_II(2, 1) = (com_I(0, 1) * com_I(0, 2)
                  - com_I(0, 0) * com_I(1, 2));

          com_II(2, 2) = (-com_I(0, 1) * com_I(0, 1) + com_I(0, 0) * com_I(1, 1));

          // get the angular velocity around the COM
          Vec com_O;

          com_O = (com_II * com_L) / denom;

          const double ekin_rot = 0.5 * com_O.dot(com_L);
          const double ekin_int = tot_ekin - ekin_trans - ekin_rot;

          trans_temp.addval(ekin_trans * 2 / (3 * gmath::physConst.get_boltzmann()));
          rot_temp.addval(ekin_rot * 2 / (3 * gmath::physConst.get_boltzmann()));

          if (constr) {
            int_temp.addval(ekin_int * 2 / ((3 * sys.mol(i).numAtoms() - 6
                    - sys.mol(i).topology().numBonds())
                    * gmath::physConst.get_boltzmann()));
          } else {
            int_temp.addval(ekin_int * 2 / ((3 * sys.mol(i).numAtoms() - 6)
                    * gmath::physConst.get_boltzmann()));
          }
        }
      }

      itrj.close();
      if (read_trv)
        itrv.close();
    }

    // print the distributions
    {
      Distribution const & a_dist = trans_temp.dist_init(100);
      cout << "# molecular translational temperature distribution\n"
              << "# min = " << trans_temp.min() << " max = " << trans_temp.max()
              << " average = " << trans_temp.ave()
              << " error estimate = " << trans_temp.ee() << "\n";
      a_dist.write_normalized(cout);
    }
    {
      Distribution const & a_dist = rot_temp.dist_init(100);
      cout << "# molecular rotational temperature distribution\n"
              << "# min = " << rot_temp.min() << " max = " << rot_temp.max()
              << " average = " << rot_temp.ave()
              << " error estimate = " << rot_temp.ee() << "\n";
      a_dist.write_normalized(cout);
    }
    {
      Distribution const & a_dist = int_temp.dist_init(100);
      cout << "# molecular internal temperature distribution\n"
              << "# min = " << int_temp.min() << " max = " << int_temp.max()
              << " average = " << int_temp.ave()
              << " error estimate = " << int_temp.ee() << "\n";
      a_dist.write_normalized(cout);
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
