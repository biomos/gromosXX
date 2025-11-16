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
#include "Hbond_calc_2c.h"

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>

#include "CubeSystem.hcc"
#include "Hbond_calc.h"
#include "groTime.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../bound/Boundary.h"

#ifdef OMP
#include <omp.h>
#endif


using utils::HB2c_calc;
using utils::HB2c;
using utils::Key2c;

void HB2c_calc::setval(gcore::System& sys, args::Arguments& args, int dummy) {
  HB_calc::setval(sys, args, dummy);
  
  Time time_args(args);
  read_time=time_args.read(); //read time from trajectory?
  time_dt=time_args.dt();
  time_start=time_args.start_time();
  
  set_reduce();
}//end HB2c_calc::setval()

void HB2c_calc::clear() {
  --frames; //-1 because in calc frames was increased by one. the ref frame is not counted
  hb2cc.clear();
  hb2cc_tmp.clear();
  ts.clear();
}//end HB2c_calc::clear()


void HB2c_calc::store_index(){
    native_key_storage.clear();
    for (HB2cContainer::const_iterator it = hb2cc_tmp.begin(); it != hb2cc_tmp.end() ; ++it) //go through native hbond map
    	native_key_storage.push_back(it->first);
}

void HB2c_calc::calc_native() {
  init_calc();
  // loop over native hydrogen bond indexes that were stored in native_key_storage
  for (size_t t = 0; t < native_key_storage.size(); ++t) {

    Key2c native_key = native_key_storage[t];
    int i,j;
    native_key.get_atom_num(i,j);

    calc(i, j);
  }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}//end HB2c_calc::calc_native()

void HB2c_calc::calc(int i, int j) {

    double dist, angles;
    gmath::Vec acceptor_j, vec_ij, bound_i;

    bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i), sys->box());
    acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j), sys->box());
    vec_ij = acceptor_j - donors.pos(i);

    if (neighbour(i, j) && distances(dist, vec_ij) && angle(i, angles, bound_i, vec_ij)) {

        Key2c key(i,j); // a unique hbond key
        hb2cc_tmp[key].add(); //for bridges calculation: use non-reduced information

        if(reduce){
            if(donors.mol(i) < 0 && !solv_donor.empty() ){ //if atom i is from a solvent molecule
                i = donors.atom(i) % donors.numSolventAtoms(); //get the solvent atom number (0,1,2,...)
                unsigned int t;
                for(t = 0; t < solv_donor.size()-1 && i != solv_donor[t]; ++t); //find the atom number of i: this is stored in solv_donor[t]. t = the position in donors
                i = solv_donor.back() + t; //solv_donor.back stores where the first solvent starts in donors: add position of the solvent in donors to where the first solvent starts in donors
            }
            if(acceptors.mol(j) < 0 && !solv_acc.empty()){ //same for acceptors
                j = acceptors.atom(j) % acceptors.numSolventAtoms();
                unsigned int t;
                for(t = 0; t < solv_acc.size()-1 && j != solv_acc[t]; ++t);
                j = solv_acc.back() + t;
            }
            key.set_index(i,j);
            hb2cc[key].add(dist, angles);
            ts.back().add_once(key);
        }
        else{
            hb2cc[key].add(dist, angles);
            ts.back().add_all(key);
        }
        ++numHb;
    }
}//end HB2c_calc::calc()

//merge hbond objects
void HB2c_calc::merge(utils::HB2c_calc& input, int traj_num){

    #ifdef OMP
    #pragma omp critical
    #endif
    {   //merge maps
        for(HB2cContainer::const_iterator it_inp = input.hb2cc.begin(); it_inp != input.hb2cc.end(); ++it_inp){
            hb2cc[it_inp->first].add(it_inp->second);
        }
        //merge other things:
        frames += input.frames; //add frames so we get the total number of frames
        traj_map[traj_num]=input.ts;
    }
}

void HB2c_calc::go_through_cubes(CubeSystem<int>& cubes_donors, CubeSystem<int>& cubes_acceptors){
    for(unsigned int c=0; c<cubes_donors.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes_donors and cubes_acceptors

        int don_atoms_i=cubes_donors.cube_i_atomlist(c).size(), //donor atoms my cube
            don_atoms_j=cubes_donors.cube_j_atomlist(c).size(), //donor atoms neighbour cube
            acc_atoms_i=cubes_acceptors.cube_i_atomlist(c).size(), //acceptor atoms my cube
            acc_atoms_j=cubes_acceptors.cube_j_atomlist(c).size(); //acceptor atoms neighbour cube

        if(don_atoms_i && acc_atoms_j) {
        //no skiP: we do not compare the same atoms in i and j!
            for(int di=0; di < don_atoms_i; ++di ){ //go through all donor atoms in first cube...
                for(int aj=0; aj < acc_atoms_j; ++aj){ //..and check against all acceptor atoms in the neighbour cube
                    calc(cubes_donors.cube_i_atomlist(c).at(di), cubes_acceptors.cube_j_atomlist(c).at(aj));
                }
            }
        }
        if(don_atoms_j && acc_atoms_i && !cubes_donors.same_cube(c)){
            for(int dj=0; dj < don_atoms_j; ++dj ){ //then do the same for the donor atoms in the second cube. this must be done since the cubes do not have redundancy. every cube pair only exists once. e.g. only once instance of 0-2 and 2-0
                for(int ai=0; ai < acc_atoms_i; ++ai){
                    calc(cubes_donors.cube_j_atomlist(c).at(dj), cubes_acceptors.cube_i_atomlist(c).at(ai));
                }
            }
        }
    }
}

void HB2c_calc::calc_hb(CubeSystem<int>& cubes_donors, CubeSystem<int>& cubes_acceptors){
    init_calc();

    if(cubes_donors.size() != cubes_acceptors.size())
        throw gromos::Exception("hbond","Donor & Acceptor CubeSystems do not have the same size. Please treat them equally!");

    //donAB & (accA + accB + accAB)
    if(!donAB.empty() && !(accA.empty() && accB.empty() && accAB.empty() ) ){
		//assign donAB atoms to cubes
		for (unsigned int i = 0; i < donAB.size(); ++i)
		    cubes_donors.assign_atom(donAB[i], donors.pos(donAB[i]));

		//assign accA + accB + accAB atoms to cubes
		for (unsigned int i = 0; i < accA.size(); ++i)
		    cubes_acceptors.assign_atom(accA[i], acceptors.pos(accA[i]));

        for (unsigned int i = 0; i < accB.size(); ++i)
		    cubes_acceptors.assign_atom(accB[i], acceptors.pos(accB[i]));

        for (unsigned int i = 0; i < accAB.size(); ++i)
		    cubes_acceptors.assign_atom(accAB[i], acceptors.pos(accAB[i]));

        go_through_cubes(cubes_donors, cubes_acceptors);
    }
//#################
    //donA & (accB + accAB)
    if(!donA.empty() && !(accB.empty() && accAB.empty() )){

        cubes_donors.delete_atoms(); //delete atoms from cubes
        cubes_acceptors.delete_atoms();

        //assign donB
        for (unsigned int i = 0; i < donA.size(); ++i)
            cubes_donors.assign_atom(donA[i], donors.pos(donA[i]));

        //assign accB atoms to cubes
        for (unsigned int i = 0; i < accB.size(); ++i)
            cubes_acceptors.assign_atom(accB[i], acceptors.pos(accB[i]));

        //assign accAB atoms to cubes
		for (unsigned int i = 0; i < accAB.size(); ++i)
		    cubes_acceptors.assign_atom(accAB[i], acceptors.pos(accAB[i]));

        go_through_cubes(cubes_donors, cubes_acceptors);
    }
//#########
    //donB & (accA + accAB)
    if(!donB.empty() && !(accA.empty() && accAB.empty() ) ){
        cubes_donors.delete_atoms(); //delete atoms from cubes
        cubes_acceptors.delete_atoms();

        //assign donB
        for (unsigned int i = 0; i < donB.size(); ++i)
            cubes_donors.assign_atom(donB[i], donors.pos(donB[i]));

        //assign accA atoms to cubes
        for (unsigned int i = 0; i < accA.size(); ++i)
            cubes_acceptors.assign_atom(accA[i], acceptors.pos(accA[i]));
        //assign accAB atoms to cubes
        for (unsigned int i = 0; i < accAB.size(); ++i)
            cubes_acceptors.assign_atom(accAB[i], acceptors.pos(accAB[i]));

        go_through_cubes(cubes_donors, cubes_acceptors);
    }

    ts.back().sort_keys(); //sort the keys in ascending order
    ts.back().set_num(numHb);
}

void HB2c_calc::calc_vac(){
    init_calc();

    //donAB & (accA + accB + accAB)
    if(!donAB.empty() && !(accA.empty() && accB.empty() && accAB.empty() ) )
        for (unsigned int i = 0; i < donAB.size(); ++i){

            for (unsigned int j = 0; j < accA.size(); ++j)
                calc(donAB[i], accA[j]);

            for (unsigned int j = 0; j < accB.size(); ++j)
                calc(donAB[i], accB[j]);

            for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donAB[i], accAB[j]);
        }

    //donA & (accB + accAB)
    if(!donA.empty() && !(accB.empty() && accAB.empty() ))
        for (unsigned int i = 0; i < donA.size(); ++i){

            for (unsigned int j = 0; j < accB.size(); ++j)
                calc(donA[i], accB[j]);

			for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donA[i], accAB[j]);
		}

	//donB & (accA + accAB)
    if(!donB.empty() && !(accA.empty() && accAB.empty() ) )
        for (unsigned int i = 0; i < donB.size(); ++i){

            for (unsigned int j = 0; j < accA.size(); ++j)
                calc(donB[i], accA[j]);

            for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donB[i], accAB[j]);
        }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}

void HB2c_calc::printstatistics(bool sort_occ, double higher){

    //open timeseries file
    opents("Hbond_2c_time_index.out", "Hbond_2c_time_numHb.out");

    cout << "# Statistics of the run:" << endl
         << endl << "# Two-centered hydrogen bonds:" << endl;

    print_header();

    unsigned int i = 1;
    for(HB2cContainer::iterator it = hb2cc.begin(); it != hb2cc.end(); ++it, ++i){
        it->second.set_id(i);
        if(it->second.num()/(double)frames * 100 >= higher)
            print(it->first);
	}

    //write timeseries - numhb and   timeseries- hbindex

    double print_time=time_start-time_dt;
    for (unsigned int traj = 0; traj < traj_map.size(); traj++) {
      for(TimeseriesContainer::const_iterator it = traj_map[traj].begin(); it != traj_map[traj].end(); ++it){
        if (!read_time) print_time = print_time+time_dt;
        else print_time = it->time();
          
        timeseriesHBtot << setw(15) << print_time
                        << setw(10) << it->num() << endl;

        for(std::vector<Key2c>::const_iterator it_key = it->keys().begin(); it_key != it->keys().end(); ++it_key){
            timeseriesHB << setw(15) << print_time
                         << setw(10) << hb2cc[*it_key].id() << endl;
        }
      }
    }

    if(sort_occ){
        HB2cContainerIteratorList hb_vec; //vector of iterators to hb2cc map

	    for(HB2cContainer::iterator it = hb2cc.begin(); it != hb2cc.end(); ++it){
    	    hb_vec.push_back(it); //populate vector with iterators to map
		}

        std::sort(hb_vec.begin(), hb_vec.end(), sort_rev_by_occ<Key2c, HB2c>);

        cout << endl << "# SORTED two-centered hydrogen bonds:" << endl;
        print_header();

		for(HB2cContainerIteratorList::const_iterator it = hb_vec.begin(); it!= hb_vec.end(); ++it){
            if((**it).second.num()/(double)frames * 100 < higher) // as soon as the occurence drops below higher value: stop
                break;
            print((**it).first);
        }

    }
}//end HB2c_calc::printstatistics()

void HB2c_calc::print_header() const{
    cout  << "#"
          << setw(6) << "HB-ID"
          << setw(7) << "Mol"
          << setw(7) << "Res"
          << setw(6) << "DONOR"
          << " -"
          << setw(7) << "Mol"
          << setw(7) << "Res"
          << setw(6) << "ACC"
          << setw(12) << "Atom"
          << setw(6) << "D"
          << " -"
          << setw(7) << "Atom"
          << setw(6) << "H"
          << " ... Atom"
          << setw(6) << "A"
          << setw(11) << "DIST"
          << setw(8) << "ANGLE"
          << setw(9) << "OCCUR"
          << setw(11) << "%" << endl;
}

void HB2c_calc::print(const Key2c& key){
    int i_d, i_a;
    key.get_atom_num(i_d,i_a);
    const HB2c& hb2cprint = hb2cc[key]; //= HB2c
    int occur = hb2cprint.num(); //how many hbonds of this index

    cout << setw(7) << hb2cprint.id() << " ";

    if (donors.mol(i_d) < 0)
        cout << setw(6) << "s"; //molecule number "s" for solvent to keep output file constistent
    else
        cout << setw(6) << donors.mol(i_d) + 1;

    cout << " "
         << setw(6) << donors.resnum(i_d) + 1 //5
         << " "
         << setw(5) << donors.resname(i_d)
         << " - ";

    if (acceptors.mol(i_a) < 0)
        cout << setw(6) << "s";
    else
        cout << setw(6) << acceptors.mol(i_a) + 1;

    cout  << " "
          << setw(6) << acceptors.resnum(i_a) + 1 //5
          << " "
          << setw(5) << acceptors.resname(i_a)
          << " "
          << setw(5+6) << bound.atom(i_d) + 1
          << " "
          << setw(5) << bound.name(i_d)
          << " - "
          << setw(6) << donors.atom(i_d) + 1
          << " "
          << setw(5) << donors.name(i_d)
          << " - "
          << setw(6) << acceptors.atom(i_a) + 1
          << " "
          << setw(5) << acceptors.name(i_a);

    cout << std::fixed
         << " "
         << setprecision(3) << setw(5+5) << hb2cprint.meandist()
         << " "
         << setprecision(3) << setw(7) << hb2cprint.meanangle()
         << " "
         << setprecision(0) << setw(8) << occur
         << " "
         << setprecision(2) << setw(10) << ((occur / (double) frames)*100)
         << endl;
}
