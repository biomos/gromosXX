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
 * @file hbond.cc
 * Monitors the occurrence of hydrogen bonds
 */

/**
 * @page programs Program Documentation
 *
 * @anchor hbond
 * @section hbond monitors the occurrence of hydrogen bonds
 * @author @ref mk, @ref ms
 * @date 9-8-2006, 13.03.2015 (ms)
 *
 * Program hbond monitors the occurrence of hydrogen bonds over molecular
 * trajectory files. It can monitor conventional hydrogen bonds, as well as
 * three-centered hydrogen bonds and solute-solvent-solute hydrogen bond bridges through geometric criteria.
 *
 * A <b>hydrogen bond</b> is considered to be present if the distance between a
 * hydrogen atom, H, connected to a donor atom D, is within a user specified
 * distance (default 0.25 nm) from an acceptor atom A and the D-H-A angle is
 * larger than the user specified angle (default 135°). If other than the default values should be used,
 * specify them with the <b>\@Hbparas</b> argument.
 *
 * Occurrences of <b>three-centered hydrogen bonds</b> are defined for a donor atom D, hydrogen
 * atom H and two acceptor atoms A1 and A2 if
 * <ol>
 * <li>the distances H-A1 and H-A2 are within a user specified value (default 0.27 nm) and </li>
 * <li>the angles D-H-A1 and D-H-A2 are larger than a second user specified value (default 90°) and </li>
 * <li> the sum of the angles D-H-A1, D-H-A2 and A1-H-A2 is larger than a third user specified value (default 340°) and</li>
 * <li> the dihedral angle defined by the planes through the atoms D-A1-A2 and H-A1-A2 is smaller than a fourth user specified value (default 15°).</li>
 * </ol>
 * To activate calculation of three-centered hydrogen bonds the flag <b>\@threecenter</b> must be given (with 4 arguments, if the default values should not be used).
 *
 * <b>Solute-solvent-solute hydrogen bond bridges</b> are calculated if the <b>\@solventbridges</b> flag is given. They report, if a solvent molecule makes a hydrogen bond
 * bridge between two solute atoms (a.k.a. protein-water-protein bridge). The geometric criteria are the same as for
 * standard hydrogen bonds.
 *
 * The user can specify two groups of atoms (A and B) between which the
 * hydrogen bonds are monitored. <b>\@DonorAtomsA</b> and <b>\@AcceptorAtomsA</b> specify the donor and acceptor atoms of group A,
 * whereas <b>\@DonorAtomsB</b> and <b>\@AcceptorAtomsB</b> specify the donor and acceptor atoms of group B. If only group A is specified, all
 * H-bonds between \@DonorAtomsA and \@AcceptorAtomsA are calculated.
 * If additionally at least one @ref AtomSpecifier of group B is given, the H-bonds are calculated as follows:
 * - \@DonorAtomsA <i>vs.</i> \@AcceptorAtomsB and
 * - \@AcceptorAtomsA <i>vs.</i> \@DonorAtomsB.
 *
 * <b>Examples:</b> (a massfile is given)
  @verbatim
 @DonorAtomsA      1:a
 @AcceptorAtomsA   1:a
 @DonorAtomsB      s:a
 @AcceptorAtomsB   s:a  @endverbatim
 * will report all H-bonds <i>between</i> molecule 1 and the solvent, whereas
  @verbatim
 @DonorAtomsA      1:a
 @AcceptorAtomsA   1:a
 @DonorAtomsB      s:a 1:a
 @AcceptorAtomsB   s:a 1:a @endverbatim
 * will additionally report all H-bonds <i>within</i> molecule 1, and
  @verbatim
 @DonorAtomsA      1:a s:a
 @AcceptorAtomsA   1:a s:a @endverbatim
 *  will report H-bonds <i>within</i> molecule 1, <i>within</i> the solvent, and <i>between</i> molecule 1 and the solvent.
 *
 * The program calculates average angles, distances and occurrences for all
 * observed hydrogen bonds over the trajectories and prints out a time series
 * of the number of observed hydrogen bonds and a time series of the H-bond IDs.<br>
 * <br><b>File names of the time series of the number of H-bonds:</b>
 * - <span style="font-family:monospace;">Hbond_2c_time_numHb.out</span> for standard H-bonds
 * - <span style="font-family:monospace;">Hbond_3c_time_numHb.out</span> for 3-centered H-bonds
 * - <span style="font-family:monospace;">Hbond_Bridges_time_numHb.out</span> for solute-solvent-solute bridge H-bonds
 *
 * <b>File names of the time series of the H-bond IDs:</b>
 * - <span style="font-family:monospace;">Hbond_2c_time_index.out</span> for standard H-bonds
 * - <span style="font-family:monospace;">Hbond_3c_time_index.out</span> for 3-centered H-bonds
 * - <span style="font-family:monospace;">Hbond_Bridges_time_index.out</span> for solute-solvent-solute bridge H-bonds
 *
 * If hydrogen bond donor and acceptor atoms
 * are not explicitly specified in the atom specifier, they can be filtered based on their masses as specified in the file given by the <b>\@massfile</b> flag.
 * The massfile must contain the blocks HYDROGENMASS and
 * ACCEPTORMASS with the masses of the hydrogens and the acceptor atoms. Optionally, a DONORMASS block can be given,
 * which contains the masses of the donor atoms (e.g. to exclude C atoms from the donors).
 *
 * If solvent is included in any AtomSpecifier, the flag <b>\@reducesolvent</b> might come in handy. It removes the (often redundant) information
 * <i>which</i> solvent molecule is part of the H-bond and thus reduces the output significantly. <i>Attention</i>: The time series of H-bond IDs
 * will not report duplicate IDs. The time series of the number of H-bonds, however, will give the total number of H-bonds.
 *
 * If a reference file is given with <b>\@ref</b>, only hydrogen bonds that are observed in the first frame of the reference file will be monitored.
 *
 * If the <b>\@sort</b> flag is given, the output will <i>additionally</i> be sorted by occurrence of the H-bonds.
 *
 * Only H-bonds that have an occurrence higher or equal than the percentage given with the <b>\@higherthan</b> flag will be printed.
 *
 * If dummy atoms are present, <b>\@excludedummy</b> flag can be used to automatically remove these atoms from the calculation.
 *
 * hbond is OMP-parallelized by trajectory file. The number of threads can be specified by the <b>\@cpus</b> flag. The number of threads will be reduced
 * to the number of trajectory files, if too many threads are requested.
 *
 * <b>Arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time"&gt; &lt;dt&gt;]  </td></tr>
 * <tr><td> \@DonorAtomsA</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@AcceptorAtomsA</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> [\@DonorAtomsB</td><td>&lt;@ref AtomSpecifier "atoms"&gt;] </td></tr>
 * <tr><td> [\@AcceptorAtomsB</td><td>&lt;@ref AtomSpecifier "atoms"&gt;] </td></tr>
 * <tr><td> [\@Hbparas</td><td>&lt;distance [nm]&gt; &lt;angle&gt; of two-centered H-bonds. Default: 0.25nm, 135°] </td></tr>
 * <tr><td> [\@threecenter</td><td>&lt;distance [nm]&gt; &lt;angle&gt; &lt;angle sum&gt; &lt;dihedral&gt; of three-centered H-bonds. Default: 0.27nm, 90.0°, 340.0°, 15.0°] </td></tr>
 * <tr><td> [\@solventbridges</td><td>Calculate solute-solvent-solute bridges] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinate file for H-bond calculation&gt;] </td></tr>
 * <tr><td> [\@massfile</td><td>&lt;massfile&gt;] </td></tr>
 * <tr><td> [\@reducesolvent</td><td>Merge output hydrogen bonds that contain solvent] </td></tr>
 * <tr><td> [\@excludedummy</td><td>If flag is used, Dummy atoms will be excluded from H-bond calculation] </td></tr>
 * <tr><td> [\@sort</td><td>Additionally print all H-bonds sorted by occurrence] </td></tr>
 * <tr><td> [\@higherthan</td><td>&lt;percentage&gt; Only print H-bonds with an occurrence higher or equal than this percentage] </td></tr>
 * <tr><td> [\@cpus</td><td>&lt;number of threads&gt; Default: 1] </td></tr>
 * <tr><td> [\@skip</td><td>skip first n frames] </td></tr>
 * <tr><td> [\@stride</td>use only every nth frame (the count is reset for every trajectory!) <td>] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * <b>Example:</b>
 * @verbatim
  hbond
    @topo             example.top
    @pbc              r
    @time             0 1
    @DonorAtomsA      1:a
    @AcceptorAtomsA   1:a
    @DonorAtomsB      s:a
    @AcceptorAtomsB   s:a
    @Hbparas          0.28 140
    @threecenter      0.30 90 330 10
    @solventbridges
    @reducesolvent
    @excludedummy
    @sort
    @higherthan       10.2
    @massfile         ../data/mass.file
    @cpus             2
    @ref              example.cnf
    @traj             example1.trc
                      example2.trc @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#ifdef OMP
#include <omp.h>
#endif

#include "../src/gromos/Exception.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Hbond_calc.h"
#include "../src/utils/Hbond_calc_2c.h"
#include "../src/utils/Hbond_calc_3c.h"
#include "../src/args/Arguments.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/System.h"
#include "../src/bound/Boundary.h"
#include "../src/bound/Triclinic.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/utils/CubeSystem.hcc"
#include "../src/utils/Hbond.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Matrix.h"


using namespace gcore;
using namespace gio;
using namespace args;
using namespace utils;
using namespace bound;
using namespace gmath;
using namespace std;

void octahedron_to_triclinic (System&, Boundary*); //for conversion from trunc oct. to triclinic

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "ref" << "DonorAtomsA" << "AcceptorAtomsA"
          << "DonorAtomsB" << "AcceptorAtomsB" << "Hbparas" << "threecenter"
          << "time" << "massfile" << "traj" << "cpus" << "gridsize" << "sort" << "solventbridges" << "reducesolvent" << "higherthan" << "skip" << "stride" << "excludedummy";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo            <molecular topology file>\n";
  usage += "\t@pbc             <boundary type>\n";
  usage += "\t[@time           <start time [ps]> <dt [ps]> Specify, if time should NOT be read from files]\n";
  usage += "\t@DonorAtomsA     <atom specifier>\n";
  usage += "\t@AcceptorAtomsA  <atom specifier>\n";
  usage += "\t[@DonorAtomsB    <atom specifier>]\n";
  usage += "\t[@AcceptorAtomsB <atom specifier>]\n";
  usage += "\t[@Hbparas        <distance [nm]> <angle [°]> Default: 0.25nm, 135°]\n";
  usage += "\t[@threecenter    <distance [nm]> <angle [°]> <minimal angle sum [°]> <dihedral [°]> Default: 0.27nm, 90.0°, 340.0°, 15.0°]\n";
  usage += "\t[@solventbridges Report solute-solvent-solute bridges]\n";
  usage += "\t[@reducesolvent  Merge output hydrogen bonds that contain solvent]\n";
  usage += "\t[@massfile       <massfile> Accepts HYDROGENMASS, ACCEPTORMASS, and DONORMASS (optional) block]\n";
  usage += "\t[@ref            <reference coordinates for H-bond calculation> Only H-bonds monitored in the first frame of this file will be reported]\n";
  usage += "\t[@sort           Additionally print all H-bonds sorted by occurrence]\n";
  usage += "\t[@excludedummy   If flag is used, Dummy atoms will be excluded from H-bond calculation]\n";
  usage += "\t[@cpus           <number of threads> Default: 1]\n";
  usage += "\t[@skip           <n> skip first n frames>]\n";
  usage += "\t[@stride         <n> use only every nth frame (per trajectory!) ]\n";
  usage += "\t[@higherthan     <percentage> Only print H-bonds with an occurrence higher or equal than this percentage]\n";
  usage += "\t@traj            <trajectory files>\n";
 // usage += "\t[@gridsize       <grid size [nm] for pairlist generation> Default: 0.5nm>]\n"; //you can use gridsize, but it will no longer be displayed for the user

  try {
#ifdef OMP
double start_total = omp_get_wtime();
double start;
#endif
    Arguments args(argc, argv, knowns, usage);
    InTopology it(args["topo"]);
    System sys(it.system());
    Time time(args);

    //required arguments that do not have a default check:
    if(args.count("pbc") <= 0){
        throw gromos::Exception("hbond", "Please specify @pbc\n\n" + usage);
    }
    if(args.count("DonorAtomsA") <= 0){
        throw gromos::Exception("hbond", "Please specify @DonorAtomsA\n\n" + usage);
    }
    if(args.count("AcceptorAtomsA") <= 0){
        throw gromos::Exception("hbond", "Please specify @AcceptorAtomsA\n\n" + usage);
    }
    if(args.count("traj") <= 0){
        throw gromos::Exception("hbond", "Please specify @traj\n\n" + usage);
    }

    Arguments::const_iterator it_arg;

    //@pbc
    string boundary_condition;
    if(args.count("pbc") > 1){
        throw gromos::Exception("hbond","@pbc: You must not give a gather method\n\n" + usage);
    }
    else{ //<= 0 was already checked
        it_arg = args.lower_bound("pbc");
        std::istringstream is(it_arg->second);
        is >> boundary_condition;
        if(boundary_condition == "t")
            cout << "# Truncated octahedral box will be converted to triclinic" << endl;
    }

    //@Hbparas
    vector<double> v_hbparas2c = args.getValues<double>("Hbparas", 2, false,
            Arguments::Default<double>() << 0.25 << 135.0);

    // get the @time argument:
    bool read_time=time.read(); //read time from trajectory?
    double time_dt=time.dt();
    double time_start=time.start_time();

    if(time_dt <= 0 || time_start < 0)
        throw gromos::Exception("hbond", "@time: Please specify <time_start> >= 0 and <dt> > 0\n\n" + usage);

    it_arg = args.lower_bound("time");

    if(args.count("time") > 2){
        std::cerr << "# Warning: hbond: @time considers only two arguments: <time_start> <dt>" << std::endl;
    }

    //load all specified atoms into atom specifier to check, if water is included: (for gridsize)
    utils::AtomSpecifier atoms(sys);

    for(it_arg=args.lower_bound("DonorAtomsA"); it_arg!=args.upper_bound("DonorAtomsA"); ++it_arg)
        atoms.addSpecifier(it_arg->second);
    if(!atoms.numSolventAtoms()) //num solvent atoms returns the number of atoms a solvent molecule has
        for(it_arg=args.lower_bound("AcceptorAtomsA"); it_arg!=args.upper_bound("AcceptorAtomsA"); ++it_arg)
            atoms.addSpecifier(it_arg->second);
    if(!atoms.numSolventAtoms())
        for(it_arg=args.lower_bound("DonorAtomsB"); it_arg!=args.upper_bound("DonorAtomsB"); ++it_arg)
            atoms.addSpecifier(it_arg->second);
    if(!atoms.numSolventAtoms())
        for(it_arg=args.lower_bound("AcceptorAtomsB"); it_arg!=args.upper_bound("AcceptorAtomsB"); ++it_arg)
            atoms.addSpecifier(it_arg->second);
    bool has_solvent = atoms.numSolventAtoms();

    //@gridsize
    double gridsize = 0.6;

    it_arg=args.lower_bound("gridsize");
    if(it_arg!=args.upper_bound("gridsize")){
        std::istringstream is(it_arg->second);
        is >> gridsize;

        if(gridsize <= v_hbparas2c[0]){
            throw gromos::Exception("hbond", "You must specify a gridsize > h-bond distance\n\n" + usage);
        }
    }

    //cout << "# Grid size: " << gridsize << " nm" << endl;

    //@solventbridges
    if(args.count("solventbridges")>=0 && !has_solvent){
        throw gromos::Exception("hbond","@solventbridges were requested, but no solvent was specified\n\n" + usage);
    }

        //@reducesolvent
    if (args.count("reducesolvent") >=0){
        if(!has_solvent)
            throw gromos::Exception("hbond","@reducesolvent was requested, but no solvent was specified\n\n" + usage);
        else
            cout << "# Pooling solvent H-bonds. NOTE: Occurences involving solvent can be greater than 100%" << endl;
    }

    //@traj
    const int traj_size = args.count("traj"); //number of trajectory files
    it_arg=args.lower_bound("traj");
    vector<Arguments::const_iterator> traj_array(traj_size);//array with pointers to trajectories:only way to go with omp

    for(int i=0; i<traj_size; ++it_arg, ++i){
        traj_array[i]=it_arg;
    }

    //@ref
    if (args.count("ref") >= 1){
        it_arg=args.lower_bound("ref");
        cout << "# Reference file: " << (it_arg->second).c_str() << endl
             << "# Notice: Only the first frame of this file will be used as reference." << endl;
    }
    if (args.count("ref") == 0){
        throw gromos::Exception("hbond","@ref: please specify a reference coordinate/trajectory file\n\n" + usage);
    }
    
    //@skip and @stride
    int stride = args.getValue<int>("stride", false, 1);
    int skip = args.getValue<int>("skip", false, 0);
    
    //@excludedummy
    //get the dummy atom IAC:
    int dummyIAC = it.forceField().dummyAtomType();
//cout << "dummy: " << dummyIAC << endl;
    if(args.count("excludedummy") >= 0)
            cout << "# Dummy atoms will be excluded" << endl;
    //check will be done in Hbond_calc.cc

    //@cpus
    int num_cpus=1;

    it_arg=args.lower_bound("cpus");
    if(it_arg!=args.upper_bound("cpus")){
        std::istringstream is(it_arg->second);
        is >> num_cpus;
        if(num_cpus <= 0)
            throw gromos::Exception("hbond","You must specify a number >0 for @cpus\n\n" + usage);
        #ifdef OMP
        if(num_cpus > traj_size){
            if(traj_size > omp_get_max_threads())
                num_cpus = omp_get_max_threads();
            else
                num_cpus = traj_size;
            cerr << "# Number of threads > number of trajectory files: not feasible. Corrected to " << num_cpus << " threads." << endl;
        }

        if(num_cpus > omp_get_max_threads()){
            cerr << "# You specified " << num_cpus << " number of threads. There are only " << omp_get_max_threads() << " threads available." << endl;
            num_cpus = omp_get_max_threads();
        }
        #else
        if(num_cpus != 1)
            throw gromos::Exception("hbond","@cpus: Your compilation does not support multiple threads. Use --enable-openmp for compilation.\n\n" + usage);
        #endif

    }
    #ifdef OMP
    omp_set_num_threads(num_cpus); //set the number of cpus for the parallel section
    #endif // OMP
    cout << "# Number of threads: " << num_cpus << endl;

    //@sort
    if (args.count("sort") >=0)
        cout << "# Requested sorting by occurrence" << endl;

    //@higherthan
    it_arg=args.lower_bound("higherthan");
    if (args.count("higherthan") > 0){
        double high;
        std::istringstream is(it_arg->second);
        is >> high;
        if(high < 0 )
            throw gromos::Exception("hbond","@higherthan argument must be >= 0\n\n" + usage);
        cout << "# Only printing H-bonds with an occurrence >= " << high << "%" << endl;
        }


    //@threecenter
    vector<double> v_hbparas3c;
    v_hbparas3c.resize(4);

    if (args.count("threecenter") >= 1) { //if some parameters were specified
      v_hbparas3c = args.getValues<double>("threecenter", 4, false,
              Arguments::Default<double>() << 0.27 << 90.0 << 340.0 << 15.0);
    }
    else if (args.count("threecenter") == 0) {
        v_hbparas3c = Arguments::Default<double>() << 0.27 << 90.0 << 340.0 << 15.0;
    }

    HBPara2c hbparas2c = HB::mk_hb2c_paras(v_hbparas2c);
    HBPara3c hbparas3c = HB::mk_hb3c_paras(v_hbparas3c);

    HB output(sys, args, hbparas2c, hbparas3c, dummyIAC); //output

    #ifdef OMP
    double prep_time = omp_get_wtime() - start_total;
    double totaltime=0;
    double calc_time_traj=0, read_time_traj=0;

    // loop over all trajectories
    #pragma omp parallel for firstprivate(sys, time) reduction(+:totaltime,calc_time_traj,read_time_traj)
    #endif
	for(int traj=0 ; traj<traj_size; ++traj){
        double frame_time = time_start - time_dt;

        HB hb(sys, args, hbparas2c, hbparas3c, dummyIAC);

        Boundary* to_pbc = new Triclinic(&sys); //in case of trunc octahedral box
        if(boundary_condition == "t")
                octahedron_to_triclinic(sys, to_pbc);

        CubeSystem<int> cubes_donors(gridsize), cubes_acceptors(gridsize);
        CubeSystem<Key2c> cubes_bridges(gridsize);

        // do native?
        if (args.count("ref") > 0) {
            InG96 ic(args["ref"]);

            if(has_solvent)
                ic.select("ALL");
            else
                ic.select("SOLUTE");

            ic >> sys;
            ic.close();

            if(boundary_condition != "v"){ //no cubesystem for vacuum
                cubes_acceptors.update_cubesystem(sys.box());
                cubes_donors = cubes_acceptors; // = cubes_acceptors; //
                cubes_bridges.update_cubesystem(sys.box());
            }
            // calculate the native hbonds:
            hb.prepare_native(cubes_donors, cubes_acceptors, cubes_bridges);
        }
        #ifdef OMP
        #pragma omp critical
        #endif
        cerr << "# Processing file: " << traj_array[traj]->second << endl;

        InG96 ic;
        // open file
        ic.open(traj_array[traj]->second);

        if(has_solvent)
            ic.select("ALL");
        else
            ic.select("SOLUTE");

        // loop over single trajectory
        int num_frames=0;
        while (!ic.eof()) {
            #ifdef OMP
            start=omp_get_wtime();
            double start_tot = omp_get_wtime();
            #endif
            if(read_time){
                ic >> sys >> time;//read coordinates & time from file
                frame_time = time.time(); //get current time
            }
            else{
                ic >> sys;
                frame_time += time_dt; //numbering starts at 0 for every traj, correct overall times are generated in printstatistics
            }
            #ifdef OMP
            read_time_traj += omp_get_wtime()-start;
            #endif

            if (!(traj==0 && num_frames < skip) && !(num_frames % stride)) {
            static int frame_check = 1;
            // get the number of atoms and break in case these numbers change from
            // one frame to another

            int numSoluAt = 0, numSolvAt = 0;
            static int numSoluAt_old = -1, numSolvAt_old = -1;
            int numSolu = sys.numMolecules();
            int numSolv = sys.numSolvents();
            for (int i = 0; i < numSolu; ++i)
              numSoluAt += sys.mol(i).numAtoms();

            for (int i = 0; i < numSolv; ++i)
              numSolvAt += sys.sol(i).numAtoms();

            if (numSoluAt_old != -1 && numSoluAt != numSoluAt_old) {
              stringstream msg;
              msg << "The number of solute atoms changed in " << traj_array[traj]->second << ":\n"
                      << "             frame " << frame_check - 1 << ": " << numSoluAt_old << " solute atoms\n"
                      << "             frame " << frame_check << ": " << numSoluAt << " solute atoms\n"
                      << "       The calculation of hbond has been stopped therefore.";
              throw gromos::Exception("hbond", msg.str());
            }
            if (numSolvAt_old != -1 && numSolvAt != numSolvAt_old) {
              stringstream msg;
              msg << "The number of solvent atoms changed in " << traj_array[traj]->second << ":\n"
                      << "             frame " << frame_check - 1 << ": " << numSolvAt_old << " solvent atoms\n"
                      << "             frame " << frame_check << ": " << numSolvAt << " solvent atoms\n"
                      << "       The calculation of hbond has been stopped therefore.";
              throw gromos::Exception("hbond", msg.str());
            }

            numSoluAt_old = numSoluAt;
            numSolvAt_old = numSolvAt;
            ++frame_check;

            hb.settime(frame_time);

            #ifdef OMP
            start=omp_get_wtime();
            #endif
            if(boundary_condition != "v"){
                cubes_acceptors.update_cubesystem(sys.box());
                cubes_donors = cubes_acceptors; // assignment operator takes over the .update_cubesystem(sys.box());
                cubes_bridges.update_cubesystem(sys.box());
            }
            #ifdef OMP
            calc_time_traj += omp_get_wtime()-start;
            start=omp_get_wtime();
            #endif

            hb.calc(cubes_donors,cubes_acceptors,cubes_bridges);
            #ifdef OMP
            calc_time_traj += omp_get_wtime()-start;
            #endif
            }
            num_frames++;
            #ifdef OMP
            totaltime += omp_get_wtime()-start_tot;
            #endif
          } // loop single traj
          
          if (traj==0 && (num_frames < skip) ) {
          cout << "#\n# WARNING: @skip is bigger than the number of frames ("<< num_frames <<")in the\n"
               << "#          first trajectory , \n"
               << "#          only the frames of the first trajectory were skipped!\n";
          }

          ic.close();
          //end of trajectory: merge hb and output. hb will be destroyed after this scope
          output.merge(hb,traj);

          delete to_pbc;
    }

    output.printstatistics();
    #ifdef OMP
    cout.precision(2);
	cout << "# Preparation time: " << prep_time << " s" << endl;
    cout << "# Total read-in time for " << traj_size << " trajectories: " << read_time_traj << " s" << endl;
    cout << "# Total calc time: " << calc_time_traj << " s" << endl;
    cout << "# Total CPU time: \t" << totaltime << " s" << endl;
    cout << "### Total real time: \t" << omp_get_wtime()-start_total << " s" << endl;
    #endif
    cout << "# hbond finished" << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void octahedron_to_triclinic (System& sys, Boundary* to_pbc){ //this bit of code is taken and modified from unify_box.cc

    Matrix rot(3,3);

    const double third = 1.0 / 3.0;
    const double sq3i = 1.0/sqrt(3.0);
    const double sq2i = 1.0/sqrt(2.0);

    const double d = 0.5*sqrt(3.0) * sys.box().K()[0];

    sys.box().K()[0] =  d;
    sys.box().K()[1] =  0.0;
    sys.box().K()[2] =  0.0;

    sys.box().L()[0] =  third * d;
    sys.box().L()[1] =  2 * third * sqrt(2.0) * d;
    sys.box().L()[2] =  0.0;

    sys.box().M()[0] = -third * d;
    sys.box().M()[1] =  third * sqrt(2.0) * d;
    sys.box().M()[2] =  third * sqrt(6.0) * d;

    sys.box().update_triclinic();

    rot = Matrix(Vec(sq3i, -2*sq2i*sq3i, 0),
	       Vec(sq3i, sq3i*sq2i, -sq2i),
	       Vec(sq3i, sq2i*sq3i, sq2i));

    gmath::Vec origo(0.0,0.0,0.0);

    for(int i=0;i<sys.numMolecules();i++){ //rotate solute
      for(int j=0;j<sys.mol(i).topology().numAtoms();j++){
        // rotate the coordinate system
        sys.mol(i).pos(j) = rot * sys.mol(i).pos(j);
        // and take the nearest image with respect to the origo
        sys.mol(i).pos(j) = to_pbc->nearestImage(origo, sys.mol(i).pos(j), sys.box());

      }
    }
    for(int j=0; j<sys.sol(0).numPos(); j++){ //rotate solvent
      // rotate the coordinate system
      sys.sol(0).pos(j) = rot * sys.sol(0).pos(j);
      // and take the nearest image with respect to the origo
      sys.sol(0).pos(j) = to_pbc->nearestImage(origo, sys.sol(0).pos(j), sys.box());
    }

    sys.box().setNtb(gcore::Box::triclinic); //reset box type
    sys.box().boxformat()=gcore::Box::genbox; //and boxformat
}
