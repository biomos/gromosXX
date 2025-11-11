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
 * @file sphericity.cc
 * Calculates the sphericity of a specified set of atoms using the moment of inertia as criterion.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sphericity
 * @section sphericity Calculates the sphericity of a specified set of atoms
 * @author @ref ms
 * @date 20-01-2018
 *
 * Calculates the sphericity of a specified set of atoms (<b>\@atoms</b>) using the principal moments of inertia (I1, I2, I3) as criterion.
 *
 * Sphericity is calculated from Tielemann et al. "Molecular Dynamics Simulations of Dodecylphosphocholine Micelles at Three Different Aggregate Sizes: Micellar Structure and Chain Relaxation" J Phys Chem B, 2000:
 * @f[ \alpha = \frac{2I_1-I_2-I_3}{I_1+I_2+I_3} @f]
 *
 * and additionally from Salaniwal et al. "Molecular Simulation of a Dichain Surfactant/Water/Carbon Dioxide System. 1. Structural Properties of Aggregates", Langmuir, 2001:
 * @f[ S=1-\frac{I_{min}}{average(I_1,I_2,I_3)} @f]
 *
 * Thereby, 0 represents a perfect sphere and 1 a highly non-spherical shape. 
 * In addition to the sphericities, the principal moments of inertia are also printed in the output.
 *
 * Sphericity is OpenMP-parallelised by trajectory files (each thread works on one <b>\@traj</b> file). Specify number of threads by <b>\@cpus</b> flag.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to calculate sphericity for&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> [\@cpus</td><td>&lt;number of CPUs to use (default: 1)&gt;]</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rot_rel
    @topo  ex.top
    @pbc   r cog
    @time  0 1
    @atoms  1:a
    @cpus 5
    @traj ex1.trc ex2.trc ex3.trc ex4.trc ex5.trc
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef OMP
#include <omp.h>
#endif

#include "../src/utils/AtomSpecifier.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Matrix.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/OutCoordinates.h"
#include "../src/gromos/Exception.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;
using namespace gmath;


int main(int argc, char **argv) {

    Argument_List knowns;
    knowns << "topo" << "pbc" << "time" << "traj" << "atoms" << "cpus"; // << "outformat";

    string usage = argv[0];
    usage += "\n\t@topo    <molecular topology file>\n";
    usage += "\t@pbc     <boundary type [gather method]>\n";
    usage += "\t@atoms     <atoms to include in sphericity calculation>\n";
    usage += "\t[@time    <time and dt>]\n";
    usage += "\t[@cpus    <numbers of CPUs to use. Default: 1>]\n";
    usage += "\t@traj    <trajectory files>\n";

    try {
        Arguments args(argc, argv, knowns, usage);

        // read topology
        InTopology it(args["topo"]);
        System sys(it.system());


        // get simulation time
        Time time(args);
        double time_dt = time.dt();
        double time_start = time.start_time();

        // parse boundary conditions



        //@traj
        const int traj_size = args.count("traj"); //number of trajectory files
        Arguments::const_iterator it_arg = args.lower_bound("traj");
        vector<Arguments::const_iterator> traj_array(traj_size);//array with pointers to trajectories:only way to go with omp

        for(int i = 0; i < traj_size; ++it_arg, ++i) {
            traj_array[i] = it_arg;
        }
        typedef map<int, vector< vector<double> > > TrajMap;
        TrajMap traj_map;

        //@cpus
        int num_cpus = 1;

        it_arg = args.lower_bound("cpus");
        if(it_arg != args.upper_bound("cpus")) {
            std::istringstream is(it_arg->second);
            is >> num_cpus;
            if(num_cpus <= 0)
                throw gromos::Exception("sphericity", "You must specify a number >0 for @cpus\n\n" + usage);
#ifdef OMP
            if(num_cpus > traj_size) {
                if(traj_size > omp_get_max_threads())
                    num_cpus = omp_get_max_threads();
                else
                    num_cpus = traj_size;
                cerr << "# Number of threads > number of trajectory files: not feasible. Corrected to " << num_cpus << " threads." << endl;
            }

            if(num_cpus > omp_get_max_threads()) {
                cerr << "# You specified " << num_cpus << " number of threads. There are only " << omp_get_max_threads() << " threads available." << endl;
                num_cpus = omp_get_max_threads();
            }
#else
            if(num_cpus != 1)
                throw gromos::Exception("sphericity", "@cpus: Your compilation does not support multiple threads. Use --enable-openmp for compilation.\n\n" + usage);
#endif

        }
#ifdef OMP
        omp_set_num_threads(num_cpus); //set the number of cpus for the parallel section
#endif // OMP
        cout << "# Number of threads: " << num_cpus << endl;

#ifdef OMP
        #pragma omp parallel for schedule (dynamic,1) firstprivate(sys, time)
#endif
        for(int traj = 0 ; traj < traj_size; ++traj) {
            double frame_time = time_start - time_dt;
            vector<double> time_vec, spher_vec, alpha_vec, I1_vec, I2_vec, I3_vec;

#ifdef OMP
            #pragma omp critical
#endif
            cout << "# Processing file: " << traj_array[traj]->second << endl;

            AtomSpecifier atoms(sys);
            for(Arguments::const_iterator
                    iter = args.lower_bound("atoms");
                    iter != args.upper_bound("atoms");
                    ++iter)
                atoms.addSpecifier(iter->second);

#ifdef OMP
            #pragma omp critical
#endif
            {
                if(atoms.empty())
                    throw gromos::Exception("sphericity", "No atoms in @atoms");
            }
            // open file
            InG96 ic;
            ic.open(traj_array[traj]->second);
            Boundary::MemPtr gathmethod; // must define variable first, otherwise compiler complains due to critical section.
            // must be critical, otherwise stdout is mangled with the couts of the gathermethod:
#ifdef OMP
            #pragma omp critical
#endif
            gathmethod = args::GatherParser::parse(sys, sys, args);
            Boundary *pbc = BoundaryParser::boundary(sys, args);
            vector<double> all_I;
            all_I.resize(3);

            // loop over single trajectory
            while(!ic.eof()) {
//                ic.select(inc);
                if(time.read()) {
                    ic >> sys >> time;//read coordinates & time from file
                    frame_time = time.time(); //get current time
                } else {
                    ic >> sys;
                    frame_time += time_dt; //numbering starts at 0 for every traj, correct overall times are generated in printstatistics
                }
                (*pbc.*gathmethod)();
                // move to com
                Vec com = PositionUtils::com(sys, atoms);
                for(int a = 0; a < atoms.size(); a++)
                    sys.mol(atoms.mol(a)).pos(atoms.atom(a)) = atoms.pos(a) - com;

                double x,y,z,xx, yy, zz, mass, spher, alpha, I1, I2, I3;
                double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;

                for(int i = 0; i < atoms.size(); ++i) {
                    Vec pos = atoms.pos(i);
                    x = pos[0]; y = pos[1]; z = pos[2];
                    xx = x*x;
                    yy = y*y;
                    zz = z*z;
                    mass = atoms.mass(i);
                    Izz += mass * (xx + yy);
                    Iyy += mass * (xx + zz);
                    Ixx += mass * (zz + yy);
                    Ixy -= mass * x * y;
                    Ixz -= mass * x * z;
                    Iyz -= mass * y * z;
                }
                Matrix inertia(Vec(Ixx,Ixy,Ixz), Vec(Ixy,Iyy,Iyz), Vec(Ixz,Iyz,Izz));
                double eigenValues[3]; bool sorting = false;
                Matrix diagonalised = inertia.diagonaliseSymmetric(eigenValues);
                I1 = eigenValues[0];
                I2 = eigenValues[1];
                I3 = eigenValues[2];

                alpha = (2*I1-I2-I3)/(I1+I2+I3);
                spher = 1 - I3 / ((I1 + I2 + I3) / 3.0);

                alpha_vec.push_back(alpha);
                spher_vec.push_back(spher);
                I1_vec.push_back(I1);
                I2_vec.push_back(I2);
                I3_vec.push_back(I3);
                time_vec.push_back(frame_time);

            }
            ic.close();
            vector< vector<double> > traj_output;
            traj_output.push_back(time_vec);
            traj_output.push_back(alpha_vec);
            traj_output.push_back(spher_vec);
            traj_output.push_back(I1_vec);
            traj_output.push_back(I2_vec);
            traj_output.push_back(I3_vec);
#ifdef OMP
            #pragma omp critical
#endif
            traj_map[traj] = traj_output;
        } // loop over traj end
        cout << "# perfect round shape: sphericity = 0. maximum deviation: sphericity=1." << endl;
        cout << "#" << setw(9) << "time" 
        << " " << setw(10) << "alpha" 
        << " " << setw(10) << "S" 
        << " " << setw(10) << "I1" 
        << " " << setw(10) << "I2" 
        << " " << setw(10) << "I3" 
        << endl;

        map<int, double > start_tme;
        for(int trj = 0; trj < traj_map.size(); ++trj) {
            if(!traj_map.count(trj))
                continue;
            start_tme[trj] = traj_map[trj][0].size() * time_dt;
        }
        double start_time = time_start;
        double tme;

        for(int trj = 0; trj < traj_map.size(); ++trj) {
            if(!traj_map.count(trj))
                continue;
            for(int n = 0; n < traj_map[trj][0].size(); ++n) {
                tme = traj_map[trj][0][n];
                if(!time.read())
                    tme += start_time;
                cout << setw(10) << tme
                     << " " << setw(10) << traj_map[trj][1][n]
                     << " " << setw(10) << traj_map[trj][2][n]
                     << " " << setw(10) << traj_map[trj][3][n]
                     << " " << setw(10) << traj_map[trj][4][n]
                     << " " << setw(10) << traj_map[trj][5][n]
                     << endl;
            }
            start_time += start_tme[trj];
        }
        cout << "# sphericity finished" << endl;
    }

    catch(const gromos::Exception &e) {
        cerr << e.what() << endl;
        exit(1);
    }
    return 0;
}
