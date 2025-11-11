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
 * @file check_box.cc
 * Check box dimensions over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor check_box
 * @section check_box Check box dimensions over a trajectory
 * @author @ref th, @ref ms
 * @date 15.07.2014
 *
 * Check_box can be used to check, if distances between atoms and periodic copies of other
 * atoms in the system get below a certain value. Check_box calculates and
 * writes out the minimum distance between any atom in the central box of
 * the system and any atom in the periodic copies (rectangular, triclinic, and truncated octahedral box are supported).
 * The gathering method must be chosen carefully, otherwise falsely low distances will be reported.
 * It is recommended to check gathering visually before using check_box.
 * If you find distances that are in the range of bond lengths, the gathering method was probably not the right one.
 *
 * Check_box is omp-parallelized by trajectory file:
 * If more threads are specified than trajectory files are given, it will correct the number of threads to match the number of trajectory files.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary condition & gather method&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;start time, dt. Specify \@time, if time should not be read from files&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to include in calculation (default: all solute - includes ions!)&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;distances below this value [nm] are reported (default: 1.4nm)&gt; </td></tr>
 * <tr><td> \@limit</td><td>&lt;number of reported distances per frame - must be less than 1000 (default: 1)&gt; </td></tr>
 * <tr><td> \@cpus</td><td>&lt;number of CPUs to use (default: 1)&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input coordinate (trajectory) files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  check_box
    @topo   example.top
    @pbc    r cog
    @time   0 0.5
    @atoms  1:1-30
    @limit  10
    @cutoff 1.3
    @cpus   2
    @traj   example1.trc
            example2.trc
 @endverbatim
 *
 * <hr>
 */

#include "../src/gromos/Exception.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#ifdef OMP
#include <omp.h>
#endif

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/CubeSystem.hcc"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Matrix.h"
#include "../src/bound/Triclinic.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

void octahedron_to_triclinic (System&, Boundary*);


struct min_dist_atoms{ //minimal distance between atoms and the involved atoms
    double mindist2;
    int min_atom_i;
    int min_atom_j;
    double time;

    //constructor
    min_dist_atoms(double md, int mai, int maj, double t):mindist2(md),min_atom_i(mai),min_atom_j(maj),time(t)
    {}
    min_dist_atoms()
    {}
};

inline bool compare_time_dist(min_dist_atoms a, min_dist_atoms b){ //sort by time and distance
    if(a.time < b.time) return true;
    if(a.time > b.time) return false;
    //time is the same: compare distances
    return a.mindist2 < b.mindist2;
}

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "traj" << "time" << "atoms" << "pbc" << "limit" << "cutoff"  << "cpus";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <periodic boundary condition> <gather method>\n";
  usage += "\t[@time      <start time [ps]> <dt [ps]> Specify, if time should NOT be read from files]\n";
  usage += "\t[@atoms     <atoms to include in calculation> Default: All solute (includes ions)]\n";
  usage += "\t[@cutoff    <atom-atom distances below this value [nm] are reported> Default: 1.4 nm>]\n";
  usage += "\t[@limit     <number of reported distances per frame> Default: 1, must be <= 1000]\n";
  usage += "\t[@cpus      <number of CPUs to use> Default: 1]\n";
  usage += "\t@traj       <input coordinate (trajectory) files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // read & check TOPO
    InTopology it(args["topo"]);

    System sys(it.system());
    System refSys(it.system());

    //init variables
    double overall_min_dist2=1e4;
    unsigned int limit = 1;
    double cutoff = 1.4;
    double cutoff2 = cutoff*cutoff;
    bool no_cutoff=false; //do not report all distances

    int num_threads=1; //number of threads used
    int num_cpus=1; //number of threads specified by the user


    //check arguments

    // get the @time argument:
    utils::Time time(args);
    bool read_time=time.read();
    double time_dt=time.dt();
    double time_start=time.start_time()-time_dt;

    //@traj
    if(args.count("traj") <= 0)
        throw gromos::Exception("check_box", "Please specify @traj");
    const int traj_size = args.count("traj"); //number of trajectory files
    Arguments::const_iterator it_arg=args.lower_bound("traj");
    vector<Arguments::const_iterator> traj_array(traj_size);//array with pointers to trajectories:only way to go with omp // now it's a vector
    for(int i=0; i<traj_size; ++it_arg, ++i){
        traj_array[i]=it_arg;
    }

    //@pbc
    if(args.count("pbc") < 2) //check that there are at least 2 arguments for pbc
            throw gromos::Exception("check_box", "@pbc: You must specify a boundary condition and a gather method");

    it_arg=args.lower_bound("pbc");
    char boundary_condition=it_arg->second.c_str()[0];

    if(boundary_condition != 'r' && boundary_condition != 'c' && boundary_condition != 't'){ //rectangular,trunc octahedron, and triclinic are supported
        stringstream msg;
        msg << "Periodic boundary of type \"" << boundary_condition << "\" not supported.";
        throw gromos::Exception("check_box", msg.str());
    }

    if(boundary_condition == 't')
        cout << "# Truncated octahedral will be converted to triclinic" << endl;


    //@limit
    it_arg=args.lower_bound("limit");
    if(it_arg!=args.upper_bound("limit")){
        std::istringstream is(it_arg->second);
        is >> limit;
        if(limit <= 0 || limit > 1000)
            throw gromos::Exception("check_box","You must specify a number >0 and <=1000 for @limit");
    }
    cout << "# limit: max. " << limit << " atom pairs per frame reported" << endl;

    //@cutoff
    it_arg=args.lower_bound("cutoff");
    if(it_arg!=args.upper_bound("cutoff")){
        std::istringstream is(it_arg->second);
        is >> cutoff;

        if(cutoff==0)
            no_cutoff=true;
        if(cutoff < 0.6)
            throw gromos::Exception("check_box","Please give a cutoff of 0 or >= 0.6");
        cutoff2 = cutoff*cutoff;
    }


    cout.precision(2);
    if(cutoff)
        cout << "# cutoff: " << cutoff << " nm" << endl;
    else{
        cout << "# cutoff: no cutoff. all distances will be reported (up to @limit)" << endl;
        cerr << "############### NOTICE ##############\n" <<
                "# You have specified a cutoff of 0: #\n"<<
                "# All distances will be calculated, #\n" <<
                "#   but this will be much slower.   #\n" <<
                "#####################################" << endl;
    }

    //@cpus
    it_arg=args.lower_bound("cpus");
    if(it_arg!=args.upper_bound("cpus")){
        std::istringstream is(it_arg->second);
        is >> num_cpus;
        if(num_cpus <= 0)
            throw gromos::Exception("check_box","You must specify a number >0 for @cpus");
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
            throw gromos::Exception("check_box","Your compilation does not support multiple threads. Use --enable-openmp for compilation.");
        #endif

    }
    #ifdef OMP
    omp_set_num_threads(num_cpus); //set the number of cpus for the parallel section
    #endif // OMP

    cout << "# number of threads: " << num_cpus << endl;

    cout.precision(2);
    typedef std::map<int, std::vector< min_dist_atoms> > TrajMap;
    TrajMap output;
    std::vector< int > num_frames;

    #ifdef OMP
    double start=omp_get_wtime();
    #pragma omp parallel for schedule (dynamic,1) firstprivate(sys,time,refSys)
    #endif
	for(int traj=0 ; traj<traj_size; ++traj){     // loop over all trajectories
        #ifdef OMP
        #pragma omp critical
        #endif
        cerr << "# Processing file: " << traj_array[traj]->second << endl; //print filename to cerr for information


        double frame_time = time_start;
        #ifdef OMP
        if(num_threads != omp_get_num_threads())
            num_threads = omp_get_num_threads();
        #endif

        std::vector<min_dist_atoms> traj_output; //each thread has it's own output vector: better parallelization

        //get the ATOMS to be included
        //this part must be here, otherwise atoms do not have coordinates
        utils::AtomSpecifier atoms(sys); //this must be initialised here, otherwise the thread's sys doesnt know the atoms

        if(args.count("atoms")>0) //if there was an atom declaration in the input
          for(Arguments::const_iterator it=args.lower_bound("atoms"); it!=args.upper_bound("atoms"); ++it)
             atoms.addSpecifier(it->second);
        else
           for(int m=0; m<sys.numMolecules(); ++m) //loop over molecules
                for(int a=0; a<sys.mol(m).numAtoms(); ++a) //loop over atoms in molecule
                    atoms.addAtom(m,a);

        int no_atoms=atoms.size();

        //a (sorted) vector to hold the minimal distance, and the involved atoms
        std::vector<min_dist_atoms> minimal_distances;
        minimal_distances.reserve(limit); //reserve memory

        std::vector<min_dist_atoms>::iterator vit;

        // Coordinates input file
        InG96 ic;
        ic.open(traj_array[traj]->second);
        ic.select("SOLUTE"); //only solute will be loaded into sys. this is default behaviour. but if something changes there...

        //get system boundaries
        Boundary::MemPtr gathmethod; // must define variable first, otherwise compiler complains due to critical section.
            // must be critical, otherwise stdout is mangled with the couts of the gathermethod:
#ifdef OMP
            #pragma omp critical
#endif
        gathmethod = args::GatherParser::parse(sys,refSys,args);
        Boundary* pbc = BoundaryParser::boundary(sys, args);
        Boundary* to_pbc = new Triclinic(&sys); //if we have a truncated octahedral, we will convert to triclinic

        CubeSystem<int> cubes(cutoff, string("OUTER"), true); //only outer neighbours, with boxshift

       // loop over all frames
        int framenum=0;
        while(!ic.eof()){
            framenum++;

            if(read_time){
                ic >> sys >> time;//read coordinates & time from file
                frame_time = time.time(); //get current time
            }
            else{
                ic >> sys;
                frame_time += time_dt; //cannot use the build in time increment, because time_start is dependend on the trajectory file number
            }

            if(boundary_condition=='t'){
                octahedron_to_triclinic(sys, to_pbc);
                (*to_pbc.*gathmethod)(); //then: gather after oct_to_tric transformation
            }
            else
                (*pbc.*gathmethod)();

            cubes.update_cubesystem(sys.box());
            //assign atoms to cubes:
            //first place the atoms so that all atoms lie right&back&up of the origin (0,0,0)
            Vec cog(0, 0, 0);
            for (int i = 0; i < no_atoms; ++i)
                cog += atoms.pos(i);
            cog /= no_atoms;
            cog -= sys.box().K()/2.0 + sys.box().L()/2.0 + sys.box().M()/2.0; //shift the cog, so that the left front lower box corner starts at (0,0,0)

            for (int i = 0; i < no_atoms; ++i){ //assign all atoms that lie within the cutoff to a cube
                atoms.pos(i) -= cog;
                cubes.assign_atom(i, atoms.pos(i)); //position all atoms around the origin of the coordinate system
            }

            for(size_t c = 0; c < cubes.size(); ++c){

                if(!cubes.cube_i_atomlist(c).empty() && !cubes.cube_j_atomlist(c).empty()){

                    Vec boxshift = cubes.boxshift(c, sys.box());
                    //if 2 neighbours have the same atoms (=there is only 1 cube total in this direction), only half of the distances must be calculated:
                    bool skip=cubes.same_cube(c);

                    for(unsigned int i=0; i < cubes.cube_i_atomlist(c).size(); ++i ){

                        int atom_i = cubes.cube_i_atomlist(c)[i];
                        Vec atom_i_pos = atoms.pos(atom_i);

                        unsigned int j=0;
                        if(skip)
                            j=i+1;

                        for(; j < cubes.cube_j_atomlist(c).size(); ++j){

                            int atom_j = cubes.cube_j_atomlist(c)[j];

                            double distance = (atom_i_pos - (atoms.pos(atom_j) + boxshift)).abs2();

                            if(distance <= cutoff2 || no_cutoff){ //if distance<=cutoff, interactions are calculated  in md++
                            //   ^^^with cutoff         ^^^no cutoff was given: all values pass
                                if(minimal_distances.size() < limit){ //fill vector
                                    for(vit = minimal_distances.begin(); vit != minimal_distances.end() && distance > vit->mindist2; ++vit); //search where to insert

                                    minimal_distances.insert(vit, min_dist_atoms(distance, atom_i, atom_j, frame_time)); //insert

                                }
                                else{ //after vector is filled: add entries and pop the last element
                                    if(distance < minimal_distances.back().mindist2){ //if the entry is smaller than the last (=biggest) element
                                        minimal_distances.pop_back(); //delete last element. this needs to be done first, otherwise the vector will be resized->slower

                                        for(vit = minimal_distances.begin(); vit != minimal_distances.end() && distance > vit->mindist2 ; ++vit); //search where to insert

                                        minimal_distances.insert(vit, min_dist_atoms(distance, atom_i, atom_j, frame_time)); //insert
                                    }
                                }
                            }
                        } //for atom j
                    } //for atom i
                }//if there are atoms in both cubes
            }//while neighbours

            if(minimal_distances.size()){ //the vector is populated
                traj_output.insert(traj_output.end(),minimal_distances.begin(),minimal_distances.end());
                    //output.push_back(minimal_distances); //this doesnt work

                //delete all list entries of this frame:
                minimal_distances.clear();
            }
        } //while frame
      num_frames.push_back(framenum);
      ic.close();

      #ifdef OMP
      #pragma omp critical
      #endif

      #ifdef OMP
      if(num_threads > 1){
        sort(traj_output.begin(), traj_output.end(), compare_time_dist);
      }
      #endif
      {
          output[traj]=traj_output;
      }
      delete to_pbc;
    } //for trajectory file //end parallel section

    utils::AtomSpecifier atoms(sys);

    if(args.count("atoms")>0) //if there was an atom declaration in the input
        for(it_arg=args.lower_bound("atoms"); it_arg!=args.upper_bound("atoms"); ++it_arg)
             atoms.addSpecifier(it_arg->second);
    else
        for(int m=0; m<sys.numMolecules(); ++m) //loop over molecules
            for(int a=0; a<sys.mol(m).numAtoms(); ++a) //loop over atoms in molecule
                atoms.addAtom(m,a);

    std::vector<min_dist_atoms>::const_iterator vec_it;
    cout.precision(2);
    cout << "###" << endl;
    cout << "# Time     Distance [nm]     Atom A - Atom B      Mol:Residue Atom A -  Mol:Resiude Atom B" << endl;
    double itime=0;
    for (unsigned int ii =0; ii < output.size(); ii++) {
      for(vec_it = output[ii].begin(); vec_it != output[ii].end(); ++vec_it){

        if (read_time) itime=vec_it->time;
        else itime=num_frames[ii] * ii * time_dt +vec_it->time;
        
        if(vec_it->mindist2 < overall_min_dist2)
            overall_min_dist2=vec_it->mindist2;

        cout << " " << left <<
                setw(12) <<
                itime <<
                setw(12) <<
                sqrt(vec_it->mindist2) << " # " <<
                right <<
                setw(7) <<
                atoms.gromosAtom(vec_it->min_atom_i)+1 << " - " <<
                left <<
                setw(7) <<
                atoms.gromosAtom(vec_it->min_atom_j)+1 <<
                " # " <<
                right <<
                setw(5) <<
                atoms.mol(vec_it->min_atom_i)+1 <<
                ":" <<
                setw(5) <<
                atoms.resnum(vec_it->min_atom_i)+1 <<
                " " <<
                left <<
                setw(5) <<
                atoms.resname(vec_it->min_atom_i) <<
                setw(4) <<
                atoms.name(vec_it->min_atom_i) <<
                "-" <<
                right <<
                setw(5) <<
                atoms.mol(vec_it->min_atom_j)+1 <<
                ":" <<
                setw(5) <<
                atoms.resnum(vec_it->min_atom_j)+1 <<
                " " <<
                left <<
                setw(5) <<
                atoms.resname(vec_it->min_atom_j) <<
                setw(4) <<
                atoms.name(vec_it->min_atom_j) <<
                endl;
      }
    }

    if(output.size()){ //if we found something
        cout << "###" << endl
             << "# OVERALL MIN: "
             << sqrt(overall_min_dist2) << " nm" << endl
             << "###" << endl;
    }
    else
        cout << "########################\n" <<
                "#  All distances were  #\n" <<
                "#   above the cutoff.  #\n" <<
                "#  Nothing to report!  #\n" <<
                "########################" << endl;
    #ifdef OMP
    cout << "# Time parallel section: " << omp_get_wtime() - start << " s" << endl;
    #endif

  }
  catch (const gromos::Exception &e){
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

    for(int i=0;i<sys.numMolecules();i++){
      for(int j=0;j<sys.mol(i).topology().numAtoms();j++){
        // rotate the coordinate system
        sys.mol(i).pos(j) = rot * sys.mol(i).pos(j);
        // and take the nearest image with respect to the origo
        sys.mol(i).pos(j) = to_pbc->nearestImage(origo, sys.mol(i).pos(j), sys.box());

      }
    }

    sys.box().setNtb(gcore::Box::triclinic);
    sys.box().boxformat()=gcore::Box::genbox;
}

