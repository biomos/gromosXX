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
 * @file diffus.cc
 * Calculates the diffusion constant for a set of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor diffus
 * @section diffus Calculates the diffusion constant 
 * @author @ref bh @ref co
 * @date 13-10-2009
 *
 * Program diffus calculates the diffusion of the
 * centre-of-geometry of a specified
 * set of atoms.
 * Firstly, the mean square displacements (@f$\Delta(t)@f$)
 * are calculated over all considered molecules and over multiple time averages.
 *
 * @f[ \Delta(t) = \frac{1}{N_m}\sum_{i=1}^{N_m}<[\vec{r_i}(t+\tau) - \vec{r_i}(\tau)]^2>_{\tau \leq t_{av}-t} @f]
 *
 * where @f$N_m@f$ is the total number of molecules (or atoms) considered in the analysis,
 * and @f$t_{av}@f$ is the duration of the averaging block.
 *
 * According to the Einstein theory, the function @f$\Delta(t)@f$ should be approximately linear and in practice,
 * the diffusion could be obtained from the slope of the @f$\Delta(t)@f$ devided by @f$2 N_d t@f$:
 *
 * @f[ D = \lim_{t\to\infty} \frac{\Delta(t)}{2 N_d t} @f]
 *
 * where @f$N_d@f$ is the number of considered dimensions (3 for 3D vectors @f$\vec{r_i}@f$).
 * The slope of the @f$\Delta(t)@f$ is obtained from linear square fit (LSF).
 * The diffus program makes an automatic LSF considering the whole time range of @f$\Delta(t)@f$, and that 
 * might not be a reasonable approach due to the poor statistics for bigger values of @f$t@f$.
 * It is strongly recommended that the user analyzes the shape of @f$\Delta(t)@f$ and perform the LSF considering 
 * only the region of linearity.
 *
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@dim</td><td>&lt;dimensions to consider&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to follow&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  diffus
    @topo  ex.top
    @pbc   r
    [@time  0 0.1]
    @dim   x y z
    @atoms s:OW
    @traj  ex.tr

   @endverbatim
 *
 * <hr>
 */
//diffus calculates diffusion
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Stat.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

// This code corresponds to an improvement of the code developed by Chris Oostenbrink (21-6-2007)
// The documentation of the old code is given below. The main improvement implemented
// in this new version is the use of a double loop for multiple time averaging, which
// provides more statistics.

// In the new code there is also no need to provide a reference structure as multiple structures
// will be considered as reference during the calculation.

// Documentation of the old code:
/* Program diffus calculates the diffusion of the centre-of-geometry of a 
 * specified set of atoms, using the Einstein equation:
 * 
 * @f[ D = \lim_{t\to\infty}\frac{<[\vec{r_0} - \vec{r}(t)]^2>}{2 N_d t} @f]
 *
 * where @f$\vec{r_0}@f$ is the centre-of-geometry in the reference
 * configuration (if none is given the program takes the first configuration of
 * the trajectory file). @f$\vec{r}(t)@f$ is the centre-of-geometry at time t.
 * @f$N_d@f$ is the number of dimensions that are being taken into account.
 *
 * The program calculates the diffusion constant by directly applying this 
 * equation as well as by applying a linear regression of the
 * mean-square-displacement. The slope of this linear regression then gives the
 * diffusion constant.
*/


int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "dim" << "atoms"  
         << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@traj   <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    // get the relevant dimensions
    int ndim = 3;
    int dim[3] = {0, 1, 2};
    {
      Arguments::const_iterator iter = args.lower_bound("dim");
      if (iter != args.upper_bound("dim")) {
        ndim = 0;

        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
        iter++;
      }
      if (iter != args.upper_bound("dim")) {
        string dum = iter->second.c_str();
        if (dum == "x") {
          dim[ndim] = 0;
          ndim++;
        } else if (dum == "y") {
          dim[ndim] = 1;
          ndim++;
        } else if (dum == "z") {
          dim[ndim] = 2;
          ndim++;
        } else {
          throw gromos::Exception("diffus",
                  "Wrong argument given for @dim !!\n"
                  "  You can give x, y, and/or z\n"
                  "  For example:\n"
                  "    @dim  x\n"
                  "  or:\n"
                  "    @dim x y z\n");
        }
      }
    }

   
    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);

    //  read reference coordinates
    System refsys(it.system());
    InG96 ic;
    int numsolv=0;


    Arguments::const_iterator iter = args.lower_bound("traj");
    if(iter != args.upper_bound("traj"))
      ic.open((iter->second).c_str());

    ic.select("ALL");
    ic >> refsys;
    numsolv = refsys.sol(0).numAtoms();
    ic.close();

    // we always need the old coordinates to take the nearest image
    System oldsys(refsys);

    // and the current system
    System sys(refsys);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refsys,args);
    // set atom number
    AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++)
        at.addSpecifier(iter->second.c_str());
    }

    // for ease of looping, we make three of these atomspecifiers
    // one for each system
    AtomSpecifier ref_at = at;
    ref_at.setSystem(refsys);
    AtomSpecifier old_at = at;
    old_at.setSystem(oldsys);


    // the reference system already contains coordinates, here we can check
    // if the ref_at actually has a size
    if (!ref_at.size())
      throw gromos::Exception("diffus",
            "No atoms to calculate the diffusion for!");

    int num_atm = at.size();
    vector<vector<Vec> > position(num_atm);
    vector<double> times;
    
    // calculate the com of the reference state
    Vec com0(0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < ref_at.size(); i++)
      com0 += ref_at.pos(i);
    com0 /= ref_at.size();

    
    int frames = 1;
    Vec comx;

    ofstream dp;
    dp.open("diffusdp.out");
    dp << "# Time series of the mean square displacement\n"
            << "# Here the calculation uses a double loop scheme:" << endl
            << "# <[r(t)-r(t+tau)]^2>/2*Nd*t\n";

    //vector<double> tdp;
    //vector<double> tt;

    frames = 0;
    // loop over all trajectories
    for(Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()) {
        ic >> sys >> time;
        if (sys.sol(0).numAtoms() != numsolv) {
           std::cerr <<  "ERROR: frame " << frames << ": number of solvents is not the same as in the reference\n";
           throw gromos::Exception("diffus", "number of solvents disagrees");
           }
        times.push_back(time.time());
        (*pbc.*gathmethod)();
        comx = Vec(0.0, 0.0, 0.0);
        Stat<double> disp_data;
        //loop over the molecules

        for(unsigned int i = 0; i < at.size(); i++) {
          at.pos(i) = pbc->nearestImage(old_at.pos(i), at.pos(i), sys.box());
          position[ i ].push_back(at.pos(i));
          old_at.pos(i) = at.pos(i);
        }
        comx /= at.size();
        frames++;
      }
      ic.close();
    }


    // This array will receive the data for final analysis

    vector<double> output_x, output_y;

    //int window=60;
    for(int it = 0; it < frames; it++) {
    //for(int it = 0; it < frames-window; it++) {
      double sum = 0;
      int counter = 0;
      double square_d = 0;
      for(int j = 0; j < frames - it; j++) {
      //for(int j = it; j < it+window; j++) {
        counter++;
        //cout << "  j : " << j << "  it : " << it << endl;
        for(unsigned int i = 0; i < at.size(); i++) {
          square_d = 0;

          for(int k = 0; k < ndim; k++) {
            const double d = (position[i][j])[dim[k]]-(position[i][j + it])[dim[k]];
            //const double d = (position[i][j])[dim[k]]-(position[i][it])[dim[k]];
            square_d += d * d;
          }

          sum += square_d;
        }
      }
      // now print out
      dp << times[it];
      //const double disp = sum / counter / at.size() / ndim;
      const double disp = sum / counter / at.size();
      dp << setw(14) << disp;
      dp << endl;
      output_y.push_back(disp);
      output_x.push_back(times[it]);
    }

    int num = output_x.size();
    double sx = 0;
    double sy = 0;
    double sxx = 0;
    double sxy = 0;


    for(int i = 0; i < num; i++) {
      sx += output_x[i];
      sy += output_y[i];
      sxx += output_x[i] * output_x[i];
      sxy += output_x[i] * output_y[i];
    }

    double a = (sxy - sx * sy / num) / (sxx - sx * sx / num);
    //double b = -(a * sx - sy) / num;

    // NOW, calculate R^2:
    double av_x = 0, av_y = 0, R = 0, nume = 0, deno1 = 0, deno2 = 0, deno3 = 0;
    av_x = sx/(output_x.size());
    av_y = sy/(output_x.size());
    for(unsigned int i=0; i<output_x.size(); i++){
      nume += (output_x[i]-av_x)*(output_y[i]-av_y);
      deno1 += pow((output_x[i]-av_x),2);
      deno2 += pow((output_y[i]-av_y),2);
    }
    deno3 = deno1*deno2;
    R = nume/pow(deno3, 0.5);
    double Rsqr = pow(R,2);
    
    cout << "# Diffusion coefficient calculated using multiple averages"<<endl;
    cout << endl;
    if (Rsqr > 0.99){
      cout << "# Diffusion was found to be linear dependent with time." << endl
              << "# This means only one regime was found and"
              << " fitting is straightforward" << endl
              << "   D = " << a / 2 / ndim 
              << "   R^2 = " << Rsqr << endl;
    }
    
    
    else{
      cout << "# ATTENTION: Number of different regimes bigger than zero or the "
              << " R^2 < 0.99" << endl
              << "# The program calculated a diffusion coefficient by linear fitting considering the whole time series of mean square displacements" << endl
              << "# However, it is STRONGLY recommended to check your diffusdp file."
              << endl
              << "   D = " << a / 2 / ndim 
              << "   R^2 = " << Rsqr << endl;

    }
    

  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
