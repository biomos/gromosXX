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
#include "Hbond_calc_bridges.h"

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

#ifdef OMP
#include <omp.h>
#endif

#include "../args/Arguments.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "Hbond_calc_2c.h"
#include "Hbond_calc.h"
#include "CubeSystem.hcc"
#include "groTime.h"



using utils::HB_bridges;
using utils::Key3c;

void HB_bridges::setval(const HB2c_calc& hb2c_calc, gcore::System& _sys, args::Arguments& _args) {
 initialise(hb2c_calc);
 sys = &_sys;
 args = &_args;
  
 Time time_args(*args);
 read_time=time_args.read(); //read time from trajectory?
 time_dt=time_args.dt();
 time_start=time_args.start_time();
 
 set_reduce();
 pbc = args::BoundaryParser::boundary(_sys, _args);
}//end HB_bridges::setval()

void HB_bridges::clear() {
  --frames; //-1 because in calc2c frames was increased by one. the ref frame is not counted
  bridges.clear();
  ts.clear();
}//end HB_bridges::clear()


//merge hbond objects
void HB_bridges::merge(const utils::HB_bridges& input, int traj_num){

    #ifdef OMP
    #pragma omp critical
    #endif
    {   //merge right map into left map
        for(BridgeContainer::const_iterator it = input.bridges.begin(); it != input.bridges.end(); ++it){ //go through rightop map
                bridges[it->first].add(it->second); //add the two ints
        }

        //merge other things:
        frames += input.frames; //add frames so we get the total number of frames
        traj_map[traj_num]=input.ts;
      }
}


void HB_bridges::calc_hb(const HB2c_calc& hb2c_calc, CubeSystem<Key2c>& cubes){
    init_calc();
    const HB2cContainer& hb2cc_tmp = hb2c_calc.get_tmp_map();

    cubes.delete_atoms();

    for (HB2cContainer::const_iterator it = hb2cc_tmp.begin(); it != hb2cc_tmp.end() ; ++it){
        int a,b;
        it->first.get_atom_num(a,b);
        Vec hb_cog = (donors.pos(a) + bound.pos(a) + acceptors.pos(b)) / 3.0; //stores the center of geometry of the 2c hbond
        cubes.assign_atom( it->first, hb_cog); //assign hbond(only the key) to a cube
    }

    for(size_t c=0; c<cubes.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes and cubes_acceptors

        int num_i=cubes.cube_i_atomlist(c).size(), //num hbonds my cube
            num_j=cubes.cube_j_atomlist(c).size(); //num hbonds neighbour cube

        if(num_i && num_j) {
            bool skip = cubes.same_cube(c); //if the cubes have the same content, we can skip half of the hbonds
            for(int i=0; i < num_i; ++i ){ //go through all hbonds in first cube...
                int j = 0;
                if(skip)
                    j=i+1;
                for(; j < num_j; ++j){ //..and check against all hbonds in the neighbour cube
                    calc(cubes.cube_i_atomlist(c).at(i), cubes.cube_j_atomlist(c).at(j));
                }
            }
        }
    }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}

void HB_bridges::calc_vac(const HB2c_calc& hb2c_calc){
    init_calc();

    const HB2cContainer& hb2cc_tmp = hb2c_calc.get_tmp_map();

    for(HB2cContainer::const_iterator it = hb2cc_tmp.begin(); it != hb2cc_tmp.end(); ++it) {
        Key2c key = it->first;
        HB2cContainer::const_iterator it_inner = it;
        if(it_inner != hb2cc_tmp.end())
            ++it_inner; //go to next 2c hb
        for(; it_inner != hb2cc_tmp.end(); ++it_inner){
            Key2c key_inner = it_inner->first;
            calc(key, key_inner);
        }
    }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}

void HB_bridges::store_index(){
    native_key_storage.clear();
    for (BridgeContainer::const_iterator it = bridges.begin(); it != bridges.end() ; ++it) //go through native hbond map
    	native_key_storage.push_back(it->first);
}

void HB_bridges::calc_native(){
    init_calc();
    for (size_t t = 0; t < native_key_storage.size(); ++t) {

        Key3c native_key = native_key_storage[t];
        int i,j,k,l;
        native_key.get_atom_num(i,j,k,l);

        calc(Key2c(i,j),Key2c(k,l));
    }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}


void HB_bridges::calc(const Key2c& key_left, const Key2c& key_right){

    int d_right, a_right, d_left, a_left; //first 2c hbond atoms

    key_left.get_atom_num(d_left,a_left); //store donor and acceptor ids
    key_right.get_atom_num(d_right,a_right);

            //only take solute - solvent - solute bridges:
    if( //A1..D1/D2..A2 and A1..D1-B1/B2-D2..A2
        (bound.atom(d_left) == bound.atom(d_right) && bound.mol(d_left) < 0 && bound.mol(d_right) < 0 && acceptors.mol(a_left) >= 0 && acceptors.mol(a_right) >= 0)  || //&& donors.atom(d_left) != donors.atom(d_right) || //only a bound atom (or maybe also the bound donor) acts as bridge
        // B1-D1..A1/A2..D2-B2
        (acceptors.atom(a_left) == acceptors.atom(a_right) && acceptors.mol(a_left) < 0 && acceptors.mol(a_right) < 0 && donors.mol(d_left) >= 0 && donors.mol(d_right) >= 0) || // an acceptor atom acts as bridge
        // B2-D2..A2/B1-D1..A1
        (bound.atom(d_left) == acceptors.atom(a_right) && bound.mol(d_left) < 0 && acceptors.mol(a_right) < 0 && acceptors.mol(a_left) >= 0 && donors.mol(d_right) >= 0) || //either outer or inner loop donor is the same atom as the inner/outer acceptor
        // B1-D1..A1/B2-D2..A2
        (bound.atom(d_right) == acceptors.atom(a_left) && bound.mol(d_right) < 0 && acceptors.mol(a_left) < 0 && acceptors.mol(a_right) >= 0 && donors.mol(d_left) >= 0)){

            Key3c key(key_left, key_right);

            if(reduce){
                if(donors.mol(d_left) < 0 && !solv_donor.empty() ){ //if atom i is from a solvent molecule
                    d_left = donors.atom(d_left) % donors.numSolventAtoms(); //get the solvent atom number (0,1,2,...)
                    unsigned int t;
                    for(t = 0; t < solv_donor.size()-1 && d_left != solv_donor[t]; ++t); //find the atom number of i: this is stored in solv_donor[t]. t = the position in donors
                    d_left = solv_donor.back() + t; //solv_donor.back stores where the first solvent starts in donors: add position of the solvent in donors to where the first solvent starts in donors
                }
                if(donors.mol(d_right) < 0 && !solv_donor.empty() ){ //if atom i is from a solvent molecule
                    d_right = donors.atom(d_right) % donors.numSolventAtoms(); //get the solvent atom number (0,1,2,...)
                    unsigned int t;
                    for(t = 0; t < solv_donor.size()-1 && d_right != solv_donor[t]; ++t); //find the atom number of i: this is stored in solv_donor[t]. t = the position in donors
                    d_right = solv_donor.back() + t; //solv_donor.back stores where the first solvent starts in donors: add position of the solvent in donors to where the first solvent starts in donors
                }
                if(acceptors.mol(a_left) < 0 && !solv_acc.empty()){ //same for acceptors
                    a_left = acceptors.atom(a_left) % acceptors.numSolventAtoms();
                    unsigned int t;
                    for(t = 0; t < solv_acc.size()-1 && a_left != solv_acc[t]; ++t);
                    a_left = solv_acc.back() + t;
                }
                if(acceptors.mol(a_right) < 0 && !solv_acc.empty()){ //same for acceptors
                    a_right = acceptors.atom(a_right) % acceptors.numSolventAtoms();
                    unsigned int t;
                    for(t = 0; t < solv_acc.size()-1 && a_right != solv_acc[t]; ++t);
                    a_right = solv_acc.back() + t;
                }
                key.set_index(d_left,a_left,d_right,a_right);
                bridges[key].add(); //create hbond and increment or just increment count by 1
                ts.back().add_once(key);
            }
            else{
                bridges[key].add(); //create hbond and increment or just increment count by 1
                ts.back().add_all(key);
            }
            ++numHb;
        }
}

void HB_bridges::printstatistics(bool sort_occ, double higher){

    //open timeseries file
    opents("Hbond_Bridges_time_index.out", "Hbond_Bridges_time_numHb.out");

    cout << endl << "# Solute-Solvent-Solute Bridges:" << endl;

    print_header();

    unsigned int i = 1;
    for(BridgeContainer::iterator it = bridges.begin(); it != bridges.end(); ++it, ++i){
        it->second.set_id(i);
        if(it->second.num()/(double)frames * 100 >= higher)
            print(it->first);
	}

    //write timeseries - numhb and   timeseries- hbindex
    double print_time=time_start-time_dt;
    for (unsigned int traj = 0; traj < traj_map.size(); traj++) {
      for(TimeseriesContainer::const_iterator it = traj_map[traj].begin(); it != traj_map[traj].end(); ++it){
        if (!read_time)  print_time = print_time+time_dt;
        else print_time = it->time();

        timeseriesHBtot << setw(15) << print_time
                        << setw(10) << it->num() << endl;

        for(std::vector<Key3c>::const_iterator it_key = it->keys().begin(); it_key != it->keys().end(); ++it_key){
            timeseriesHB << setw(15) << print_time
                         << setw(10) << bridges[*it_key].id() << endl;
        }
      }
    }

    if(sort_occ){
        BridgeContainerIteratorList hb_vec; //vector of iterators to hb2cc map

        for(BridgeContainer::iterator it = bridges.begin(); it != bridges.end(); ++it){
    	    hb_vec.push_back(it); //populate vector with iterators to map
		}

        std::sort(hb_vec.begin(), hb_vec.end(), sort_rev_by_occ<Key3c, Bridge>);

        cout << endl << "# SORTED Solute-Solvent-Solute bridges:" << endl;
        print_header();

		for(BridgeContainerIteratorList::const_iterator it = hb_vec.begin(); it!= hb_vec.end(); ++it){
            if((**it).second.num()/(double)frames * 100 < higher) // as soon as the occurence drops below higher value: stop
                break;
            print((**it).first);
        }
    }

}//end HB_bridges::printstatistics()

void HB_bridges::print_header() const {
    cout  << "#"
          << right
          << setw(7) << "HB"
          << setw(21) << "Mol:Atom Solute 1"
          << " ... "
          << setw(18) << "Bridging Solvent  "
          << " ... "
          << left
          << setw(20) << "Solute 2"
          << " # "
          << right
          << setw(20) << "ResNum Name Solute 1"
          << " ... "
          << setw(20) << "Bridging Solvent  "
          << " ... "
          << left
          << setw(20) << "Solute 2"
          << right
          << setw(10) << "OCCUR"
          << setw(11) << "%" << endl;
}

void HB_bridges::print(const Key3c& key){

    int d1, a1, d2, a2;
    key.get_atom_num(d1, a1, d2, a2);
    const Bridge& bridgeprint = bridges[key];
    int occur = bridgeprint.num(); //how many hbonds of this index

    std::ostringstream atom_num1, atom_num2, atom_num3, atom_name1, atom_name2, atom_name3;

    //first case A..HWx..A
    if(bound.atom(d1) == bound.atom(d2) && donors.atom(d1) == donors.atom(d2)){
        atom_num1 << acceptors.mol(a1)+1 << ':' << acceptors.atom(a1)+1;
        atom_num2 << "S:" << donors.atom(d1)+1;
        atom_num3 << acceptors.mol(a2)+1 << ':' << acceptors.atom(a2)+1;

        atom_name1 << acceptors.resnum(a1)+1 << ' ' << acceptors.resname(a1) << ' ' << acceptors.name(a1);
        atom_name2 << donors.resnum(d1)+1 << ' ' << donors.name(d1);
        atom_name3 << acceptors.resnum(a2)+1 << ' ' << acceptors.resname(a2) << ' ' << acceptors.name(a2);
    }
    //second case: D-H..OW..H-D
    else if(acceptors.atom(a1) == acceptors.atom(a2)){
        atom_num1 << donors.mol(d1)+1 << ':' << bound.atom(d1)+1 << '-' << donors.atom(d1)+1;
        atom_num2 << "S:" << acceptors.atom(a1)+1;
        atom_num3 << donors.mol(d2)+1 << ':' << donors.atom(d2)+1 << '-' << bound.atom(d2)+1;

        atom_name1 << donors.resnum(d1)+1 << ' ' << donors.resname(d1) << ' ' << bound.name(d1) << '-' << donors.name(d1);
        atom_name2 << acceptors.resnum(a1)+1 << ' ' << acceptors.name(a1);
        atom_name3  << donors.resnum(d2)+1 << ' ' << donors.resname(d2) << ' ' << donors.name(d2) << '-' << bound.name(d2);
    }
    //third case: A..HW1-OW-HW2..A
    else if(bound.atom(d1) == bound.atom(d2) && donors.atom(d1) != donors.atom(d2)){
        atom_num1 << acceptors.mol(a1)+1 << ':' << acceptors.atom(a1)+1;
        atom_num2 << "S:" << donors.atom(d1)+1 << '-' << bound.atom(d2)+1 << '-' << donors.atom(d2)+1;
        atom_num3 << acceptors.mol(a2)+1 << ':' << acceptors.atom(a2)+1;

        atom_name1 << acceptors.resnum(a1)+1 << ' ' << acceptors.resname(a1) << ' ' << acceptors.name(a1);
        atom_name2 << donors.resnum(d1)+1 << ' ' << donors.name(d1) << '-' << bound.name(d1) << '-' << donors.name(d2);
        atom_name3 << acceptors.resnum(a2)+1 << ' ' << acceptors.resname(a2) << ' ' << acceptors.name(a2);
    }
    //fourth case: D-H..OW-HWx..A
    else if(acceptors.atom(a1) == bound.atom(d2) || acceptors.atom(a2) == bound.atom(d1)){
        if(acceptors.mol(a1) < 0){
            atom_num1 << donors.mol(d1)+1 << ':' << bound.atom(d1)+1 << '-' << donors.atom(d1)+1;
            atom_num2 << "S:" << bound.atom(d2)+1 << '-' << donors.atom(d2)+1;
            atom_num3 << acceptors.mol(a2)+1 << ':' << acceptors.atom(a2)+1;
            atom_name1 << donors.resnum(d1)+1 << ' ' << donors.resname(d1) << ' ' << bound.name(d1) << '-' << donors.name(d1);
            atom_name2 << donors.resnum(d2)+1 << ' ' << bound.name(d2) << '-' << donors.name(d2);
            atom_name3 << acceptors.resnum(a2)+1 << ' ' << acceptors.resname(a2) << ' ' << acceptors.name(a2);
        }
        else if(acceptors.mol(a2) < 0){
            atom_num1 << donors.mol(d2)+1 << ':' << bound.atom(d2)+1 << '-' << donors.atom(d2)+1;
            atom_num2 << "S:" << bound.atom(d1)+1 << '-' << donors.atom(d1)+1;
            atom_num3 << acceptors.mol(a1)+1 << ':' << acceptors.atom(a1)+1;
            atom_name1 << donors.resnum(d2)+1 << ' ' << donors.resname(d2) << ' ' << bound.name(d2) << '-' << donors.name(d2);
            atom_name2 << donors.resnum(d1)+1 << ' ' << bound.name(d1) << '-' << donors.name(d1);
            atom_name3 << acceptors.resnum(a1)+1 << ' ' << acceptors.resname(a1) << ' ' << acceptors.name(a1);
        }
    }

    int pad=20, padding1,padding2;
    padding1 = (pad - atom_num2.str().length())/2;
    padding2 = pad - padding1 - atom_num2.str().length();
    cout      << setw(8) << bridgeprint.id()
              << " "
              << right
              << setw(pad) << atom_num1.str()
              << " ..."
              << setw(padding1) << " "
              << atom_num2.str()
              << setw(padding2) << " "
              << "... "
              << left
              << setw(pad) << atom_num3.str()
              << " # ";
    padding1 = (pad - atom_name2.str().length())/2;
    padding2 = pad - padding1 - atom_name2.str().length();
    cout << right
              << setw(pad) << atom_name1.str()
              << " ... "
              << setw(padding1) << " "
              << atom_name2.str()
              << setw(padding2) << " "
              << " ... "
              << left
              << setw(pad) << atom_name3.str();

    cout.precision(0);
    cout << right << setw(10) << occur << " ";
    cout.setf(ios::floatfield, ios::fixed);
    cout.precision(2);
    cout << setw(10) << ((occur / (double) frames)*100)
         << endl;
}
