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
#include "Hbond_calc_3c.h"

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>

#include "../args/Arguments.h"
#include "../bound/Boundary.h"
#include "../args/BoundaryParser.h"
#include "../gromos/Exception.h"
#include "CubeSystem.hcc"
#include "Hbond_calc.h"
#include "Hbond_calc_2c.h"
#include "groTime.h"

using utils::HB3c_calc;
using utils::HB3c;
using utils::Key3c;


void HB3c_calc::setval(const HB2c_calc& hb2c_calc, gcore::System& _sys, args::Arguments& _args) {
  sys = &_sys;
  args = &_args;
  initialise(hb2c_calc);
  
  Time time_args(*args);
  read_time=time_args.read(); //read time from trajectory?
  time_dt=time_args.dt();
  time_start=time_args.start_time();
  
  set_reduce();
  pbc = args::BoundaryParser::boundary(_sys, _args);
}//end HB3c_calc::setval()

void HB3c_calc::clear() {
  --frames;
  hb3cc.clear();
  hb3cc_tmp.clear();
  ts.clear();
}//end HB3c_calc::clear()

void HB3c_calc::store_index(){
    native_key_storage.clear();
    for (HB3cContainer::const_iterator it = hb3cc_tmp.begin(); it != hb3cc_tmp.end() ; ++it) //go through native hbond map
    	native_key_storage.push_back(it->first);
}

void HB3c_calc::calc_native() {
  init_calc();

  for (size_t t = 0; t < native_key_storage.size(); ++t) {
    Key3c native_index = native_key_storage[t];

    int i,j,k,l;
    native_index.get_atom_num(i,j,l,k);
    std::vector<int> k_vec(1, k); //k must be a vector. in this case only with 1 element

    calc(i, j, k_vec);
  }
  ts.back().sort_keys();
  ts.back().set_num(numHb);
}//end HB3c_calc::calc_native()

void HB3c_calc::calc(int i, int j, const std::vector<int>& k_atoms) {

  double angle1, dist1;
  gmath::Vec acceptor_j, bound_i, vec_ij;

  bound_i = pbc->nearestImage(donors.pos(i), bound.pos(i), sys->box());
  acceptor_j = pbc->nearestImage(donors.pos(i), acceptors.pos(j), sys->box());
  vec_ij = acceptor_j - donors.pos(i);

  if (neighbour(i, j) && distances(dist1, vec_ij) && angle(i, angle1, bound_i, vec_ij)) {

	//std::vector<int>::const_iterator k_it = std::find(k_atoms.begin(), k_atoms.end(), j);
	//if(k_it == k_atoms.end()) //j was not found in atoms_k
		//k_it = k_atoms.begin();

    for(std::vector<int>::const_iterator k_it = k_atoms.begin(); k_it != k_atoms.end(); ++k_it){

      int k = *k_it;

      double angle2, dist2;
      gmath::Vec acceptor_k, vec_ik;

      acceptor_k = pbc->nearestImage(donors.pos(i), acceptors.pos(k), sys->box());
      vec_ik = acceptor_k - donors.pos(i);

      //always create the key such that a1 < a2; this way we only need to search hb3cc_tmp (current frame only),
      //and not hb3cc (all frames, much bigger!!), to check if the key already exists.
      //this way we dont need to search for the inverse key, which must be searched in hb3cc (slow!)
      Key3c key;
      if(j<=k)
        key.set_index(i,j,i,k); //d1,a1,d2,a2; d1==d2
      else
        key.set_index(i,k,i,j);

		//             vvvv if we found this hbond in this frame already
      if (j != k && !hb3cc_tmp.count(key) && neighbour(i, k) && distances(dist2, vec_ik) && angle(i, angle2, bound_i, vec_ik)) {

        double angle_sum, dihedral;
        angle_sum = angle1 + angle2;

        if( anglesum(i, angle_sum, acceptor_j, acceptor_k) &&
            dihedrals(i, dihedral, bound_i, acceptor_j, acceptor_k)) {

            hb3cc_tmp[key].add(); //add to tmp map. each hbond can only occur once per frame

            if(reduce){

                if(donors.mol(i) < 0){ //if atom i is from a solvent molecule
                    i = donors.atom(i) % donors.numSolventAtoms(); //get the solvent atom number (0,1,2,...)
                    unsigned int t;
                    for(t = 0; t < solv_donor.size()-1 && i != solv_donor[t]; ++t); //find the atom number of i: this is stored in solv_donor[t]. t = the position in donors
                    i = solv_donor.back() + t; //solv_donor.back stores where the first solvent starts in donors: add position of the solvent in donors to where the first solvent starts in donors
                }
                if(acceptors.mol(j) < 0){ //same for acceptors
                    j = acceptors.atom(j) % acceptors.numSolventAtoms();
                    unsigned int t;
                    for(t = 0; t < solv_acc.size()-1 && j != solv_acc[t]; ++t);
                    j = solv_acc.back() + t;
                }
                if(acceptors.mol(k) < 0){
                    k = acceptors.atom(k) % acceptors.numSolventAtoms();
                    unsigned int t;
                    for(t = 0; t < solv_acc.size()-1 && k != solv_acc[t]; ++t);
                    k = solv_acc.back() + t;
                }
                key.set_index(i,j,i,k); //i,j,k could have been changed ^^^^

                hb3cc[key].add(dist1, dist2, angle1, angle2, angle_sum, dihedral); //add to complete map
                ts.back().add_once(key);
            }
            else{
                hb3cc[key].add(dist1, dist2, angle1, angle2, angle_sum, dihedral);
                ts.back().add_all(key);
            }

            ++numHb;
          }
      }
    }
  }
}//end HB3c_calc::calc()

void HB3c_calc::go_through_cubes(CubeSystem<int>& cubes_donors, CubeSystem<int>& cubes_acceptors){

    for(unsigned int c=0; c<cubes_donors.size(); ++c){ //go through all cubes with donor atoms from A and all cubes with acceptor atoms from B. the size() is the same for cubes_donors and cubes_acceptors

        int don_atoms_i=cubes_donors.cube_i_atomlist(c).size(), //donor atoms my cube
            don_atoms_j=cubes_donors.cube_j_atomlist(c).size(), //donor atoms neighbour cube
            acc_atoms_i=cubes_acceptors.cube_i_atomlist(c).size(), //acceptor atoms my cube
            acc_atoms_j=cubes_acceptors.cube_j_atomlist(c).size(); //acceptor atoms neighbour cube

        std::vector<int> k_list;
        if((don_atoms_i && acc_atoms_j) || (don_atoms_j && acc_atoms_i)) //only update if needed
            k_list = cubes_acceptors.neighbour_atomlist(c); //cubes_k_list atoms k list

        if(don_atoms_i && acc_atoms_j) {
        //no skiP: we do not compare the same atoms in i and j!
            for(int di=0; di < don_atoms_i; ++di ){ //go through all donor atoms in first cube...
                for(int aj=0; aj < acc_atoms_j; ++aj){ //..and check against all acceptor atoms in the neighbour cube
                    calc(cubes_donors.cube_i_atomlist(c).at(di), cubes_acceptors.cube_j_atomlist(c).at(aj), k_list);
                }
            }
        }
        if(don_atoms_j && acc_atoms_i && !cubes_donors.same_cube(c)){
            for(int dj=0; dj < don_atoms_j; ++dj ){ //then do the same for the donor atoms in the second cube. this must be done since the cubes do not have redundancy. every cube pair only exists once. e.g. only once instance of 0-2 and 2-0
                for(int ai=0; ai < acc_atoms_i; ++ai){
                    calc(cubes_donors.cube_j_atomlist(c).at(dj), cubes_acceptors.cube_i_atomlist(c).at(ai), k_list);
                }
            }
        }
        //no need to clear k_list, it is refilled every time it is needed
    }
}

void HB3c_calc::calc_hb(CubeSystem<int>& cubes_donors, CubeSystem<int>& cubes_acceptors){
    init_calc();
//#################

    if(cubes_donors.size() != cubes_acceptors.size())
        throw gromos::Exception("hbond","Donor & Acceptor CubeSystems do not have the same size. Please treat them equally!");


    //donAB & (accA + accB + accAB)
    if(!donAB.empty() && !(accA.empty() && accB.empty() && accAB.empty() ) ){

        cubes_donors.delete_atoms(); //delete atoms from cubes
        cubes_acceptors.delete_atoms();

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


    ts.back().sort_keys();
    ts.back().set_num(numHb);
}

void HB3c_calc::calc_vac(){
    init_calc();

    std::vector<int> k_list;

    //donAB & (accA + accB + accAB)
    if(!donAB.empty() && !(accA.empty() && accB.empty() && accAB.empty() ) ){
        k_list = accA;
        k_list.insert(k_list.end(), accB.begin(), accB.end());
        k_list.insert(k_list.end(), accAB.begin(), accAB.end());

        for (unsigned int i = 0; i < donAB.size(); ++i){

            for (unsigned int j = 0; j < accA.size(); ++j)
                calc(donAB[i], accA[j], k_list);

            for (unsigned int j = 0; j < accB.size(); ++j)
                calc(donAB[i], accB[j], k_list);

            for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donAB[i], accAB[j], k_list);
        }
    }


    //donA & (accB + accAB)
    if(!donA.empty() && !(accB.empty() && accAB.empty() )){
        k_list = accB; //assignment deletes everything that was in the vector before
        k_list.insert(k_list.end(), accAB.begin(), accAB.end());
        for (unsigned int i = 0; i < donA.size(); ++i){

            for (unsigned int j = 0; j < accB.size(); ++j)
                calc(donA[i], accB[j], k_list);

			for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donA[i], accAB[j], k_list);
		}
    }

	//donB & (accA + accAB)
    if(!donB.empty() && !(accA.empty() && accAB.empty() ) ){
        k_list = accA; //assignment deletes everything that was in the vector before
        k_list.insert(k_list.end(), accAB.begin(), accAB.end());
        for (unsigned int i = 0; i < donB.size(); ++i){

            for (unsigned int j = 0; j < accA.size(); ++j)
                calc(donB[i], accA[j], k_list);

            for (unsigned int j = 0; j < accAB.size(); ++j)
                calc(donB[i], accAB[j], k_list);
        }
    }
    ts.back().sort_keys();
    ts.back().set_num(numHb);
}

void HB3c_calc::printstatistics(bool sort_occ, double higher){

    opents("Hbond_3c_time_index.out", "Hbond_3c_time_numHb.out");
    cout << endl << "# Three-centered hydrogen bonds:" << endl;

    print_header();

    unsigned int i = 1;
    for(HB3cContainer::iterator it = hb3cc.begin(); it != hb3cc.end(); ++it, ++i){
        it->second.set_id(i);
        if(it->second.num()/(double)frames * 100 >= higher)
            print(it->first);
	}

    //write timeseries - numhb and   timeseries- hbindex
    double print_time=time_start-time_dt;
    for (unsigned int traj = 0; traj < traj_map.size(); traj++) {
      for(TimeseriesContainer::const_iterator it = traj_map[traj].begin(); it != traj_map[traj].end(); ++it){
        if (!read_time) {
          print_time = print_time+time_dt;
        }
        else print_time = it->time();

        timeseriesHBtot << setw(15) << print_time
                        << setw(10) << it->num() << endl;

        for(std::vector<Key3c>::const_iterator it_key = it->keys().begin(); it_key != it->keys().end(); ++it_key){
            timeseriesHB << setw(15) << print_time
                         << setw(10) << hb3cc[*it_key].id() << endl;
        }
      }
    }

    if(sort_occ){
        HB3cContainerIteratorList hb_vec; //vector of iterators to hb2cc map

	    for(HB3cContainer::iterator it = hb3cc.begin(); it != hb3cc.end(); ++it){
    	    hb_vec.push_back(it); //populate vector with iterators to map
		}

        std::sort(hb_vec.begin(), hb_vec.end(), sort_rev_by_occ<Key3c, HB3c>);

        cout << endl << "# SORTED three-centered hydrogen bonds:" << endl;
        print_header();

		for(HB3cContainerIteratorList::const_iterator it = hb_vec.begin(); it!= hb_vec.end(); ++it){
            if((**it).second.num()/(double)frames * 100 < higher) // as soon as the occurence drops below higher value: stop
                break;
            print((**it).first);
        }
    }

}//end HB3c_calc::printstatistics()

void HB3c_calc::print_header() const{
    cout  << right
          << "#"
          << setw(6) << "HB-ID"
          << setw(7) << "Mol"
          << setw(7) << "Res"
          << setw(6) << "DONOR"
          << " -"
          << setw(7) << "Mol"
          << setw(7) << "Res"
          << setw(6) << "ACC1"
          << setw(12) << "Atom"
          << setw(6) << "D"
          << " -"
          << setw(7) << "Atom"
          << setw(6) << "H"
          << " ... Atom"
          << setw(6) << "A1"
          << setw(11) << "DIST1"
          << setw(8) << "ANGLE1"
          << setw(9) << "ANGLESUM"
          << setw(9) << "DIHED."
          << setw(9) << "OCCUR"
          << setw(11) << "%"
          << endl;
          //second line:

    cout << "#"
          << setw(28) << "\\"
          << setw(7) << "Mol"
          << setw(7) << "Res"
          << setw(6) << "ACC2"
          << setw(34) << " \\"
          << "... Atom"
          << setw(6) << "A2"
          << setw(11) << "DIST2"
          << setw(8) << "ANGLE2"
          << endl;

}

void HB3c_calc::print(const Key3c& key){

    int i_d, i_a1, i_a2,i_d2;
    key.get_atom_num(i_d,i_a1,i_d2,i_a2);
    const HB3c& hb3cprint = hb3cc[key]; //= HB3c

    cout << setw(7) << hb3cc[key].id() << " ";

    if (donors.mol(i_d) < 0)
        cout << setw(6) << "s";
    else
        cout << setw(6) << donors.mol(i_d) + 1;

    cout << " "
         << setw(6) << donors.resnum(i_d) + 1 //5
         << " "
         << setw(5) << donors.resname(i_d)
         << " - ";

    if (acceptors.mol(i_a1) < 0)
        cout << setw(6) << "s";
    else
        cout << setw(6) << acceptors.mol(i_a1) + 1;

    cout  << " "
          << setw(6) << acceptors.resnum(i_a1) + 1 //5
          << " "
          << setw(5) << acceptors.resname(i_a1)
          << " "
          << setw(5+6) << bound.atom(i_d) + 1
          << " "
          << setw(5) << bound.name(i_d)
          << " - "
          << setw(6) << donors.atom(i_d) + 1
          << " "
          << setw(5) << donors.name(i_d)
          << " - "
          << setw(6) << acceptors.atom(i_a1) + 1
          << " "
          << setw(5) << acceptors.name(i_a1);

    cout << std::fixed
        << setprecision(3)
        << " "
        << setw(5+5) << hb3cprint.meandist(0)
        << " "
        << setw(7) << hb3cprint.meanangle(0)
        << " "
        << setw(8) << hb3cprint.meanangle_sum()
        << " "
        << setw(8) << hb3cprint.meandihedral()
        << " "
        << setprecision(0) << setw(8) << hb3cprint.num()
        << " "
        << setprecision(2) << setw(10) << ((hb3cprint.num() / (double) frames)*100)
        << endl;

    // and the second line
    cout << setw(30) << "\\ ";

    if (acceptors.mol(i_a2) < 0)
        cout << setw(6) << " ";
    else
        cout << setw(6) << acceptors.mol(i_a2) + 1;

    cout  << " "
          << setw(6) << acceptors.resnum(i_a2) + 1 //5
          << " "
          << setw(5) << acceptors.resname(i_a2)
          << setw(36) << "\\ "
          << setw(6) << acceptors.atom(i_a2) + 1
          << " "
          << setw(5) << acceptors.name(i_a2)
          << " "
          << setprecision(3)
          << setw(10) << hb3cprint.meandist(1)
          << " "
          << setw(7) << hb3cprint.meanangle(1)
          << endl;
}

//merge hbond objects
void HB3c_calc::merge(utils::HB3c_calc& input, int traj_num){
    #ifdef OMP
    #pragma omp critical
    #endif
    {   //merge maps:
        for(HB3cContainer::const_iterator it=input.hb3cc.begin(); it != input.hb3cc.end(); ++it){ //go through input map
                hb3cc[it->first].add(it->second); //add the HB3c entries
        }

        //merge other things:
        frames += input.frames; //add frames so we get the total number of frames
        traj_map[traj_num]=input.ts;
      }
}

