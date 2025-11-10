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
 * @file eps_field.cc
 * Calculate the relative permittivity over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor eps_field
 * @section eps_field Calculate the relative permittivity when an external electric field was used
 * @author @ref sr
 * @date 7.9.10
 *
 * Program eps_field estimates the static dielectric permittivity,
 * @f$\epsilon(0)@f$, of a liquid when an external electric field was applied
 * during the simulation. The permittivity for a specific external field is given by
 *
 * @f[ \epsilon(0) = 1 + 4 \pi \frac{<P_{z}>_{t}}{E^{ex}_{z}} @f]
 *
 * where @f$E^{ex}@f$ is the external electric field,
 * @f$\epsilon_0@f$ is the dielectric permittivity of vacuum,
 * and @f$\vec{P}@f$ is the polarisation of the system defined as
 *
 * @f[ \vec{P}(t) = V(t)^{-1}  \vec{M}(t) @f]
 *
 * where @f$\vec{M}@f$ is the total dipole moment of the system and @f$V@f$ is
 * the volume.
 *
 * Note to get a linear response of the polarisation, the electric field should
 * be small enough to avoid saturation, which is the case if
 *
 * @f[ \frac{<\mu_{i,z} E^{ex}_{z}>}{3k_{B}} << T @f]
 *
 * with @f$\mu_{i}@f$ the dipole moment of molecule @f$i@f$,
 * @f$k_{B}@f$ the Boltzmann constant and @f$T@f$ the temperature, is fulfilled.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@E_ex</td><td>&lt;z-component of external electric field&gt; </td></tr>
 * <tr><td> \@trs</td><td>&lt;special trajectory files (with box dipole moments)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  eps_field
    @topo  ex.top
    @pbc   r
    [@time  0 0.2]
    @E_ex  0.1
    @trs   ex.trs
 @endverbatim
 *
 * <hr>
 */

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <cstdlib>

#include <math.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/utils/DipTraj.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "E_ex" << "trs";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pbc     <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@E_ex    <external electric field in z-direction>\n";
  usage += "\t@trs     <special trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // read E_ex
    double E_ex = args.getValue<double>("E_ex");
    cout << "#Note: Only the dipole moment in z-direction was considered for calculation." << endl;

    // parse boundary conditions
    //Boundary *pbc = BoundaryParser::boundary(sys, args);

    // prepare the calculation of the average volume
    //double vcorr = 1.0;
    //if (pbc->type() == 't')
    //  vcorr = 0.5;

    // and values to calculate the dipole fluctuation
    double sum_pz = 0.0, sum_px = 0.0, sum_py = 0.0;
    int numFrames = 0;

    // define input trajectory
    DipTraj boxdip;

    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("trs"),
      to = args.upper_bound("trs"); iter != to; ++iter) {

      // open file
      boxdip.open((iter->second).c_str());

      // loop over single trajectory
      while (!boxdip.eof()) {
        boxdip.read();

        DipData dip;
        boxdip >> dip >> time;
        
        gmath::Vec p(0.0,0.0,0.0);
        for (unsigned int i = 0; i < dip.data().size(); ++i) {
          p = dip.data()[i].dipole / dip.data()[i].vol;          
          sum_px += p[0];
          sum_py += p[1];
          sum_pz += p[2];
        }

        cout << setw(15) << time.time() << setw(15) << p[2] << endl;

        numFrames++;
      }
      boxdip.close();
    }
    sum_pz /= numFrames;
    sum_px /= numFrames;
    sum_py /= numFrames;

    double eps0 = 1 + 4*gmath::physConst.get_pi() *(sum_pz/E_ex);

    cout.precision(8);
    cout << "# " << setw(15) << "Px" << setw(15) << "Py" << setw(15) << "Pz"
            << setw(15) << "e(0)" << endl;
    cout << "# " << setw(15) << sum_px << setw(15) << sum_py << setw(15) << sum_pz
            << setw(15) << eps0 << endl;
    

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

