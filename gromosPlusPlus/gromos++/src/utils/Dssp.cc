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
#include "Dssp.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <vector>

#include "AtomSpecifier.h"
#include "Neighbours.h"
#include "../args/Arguments.h"
#include "../args/BoundaryParser.h"
#include "../args/GatherParser.h"
#include "../bound/Boundary.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gmath/Vec.h"

using namespace args;
using namespace gio;
using namespace gcore;
using namespace bound;
using namespace std;

using gcore::System;
using args::Arguments;
using utils::Dssp;
using utils::AtomSpecifier;

void Dssp::determineAtoms(utils::AtomSpecifier &protein) {
  protein.sort();
  for (unsigned int m = 1; m < protein.size(); m++) {
    if (protein.mol(m) == protein.mol(m - 1)) {
      if (protein.name(m - 1) == "N" && protein.name(m) == "H") {
	d_H.addAtom(protein.mol(m), protein.atom(m));
	d_N.addAtom(protein.mol(m - 1), protein.atom(m - 1));
      }
      if (protein.name(m - 1) == "C" && protein.name(m) == "O") {
	d_O.addAtom(protein.mol(m), protein.atom(m));
	d_C.addAtom(protein.mol(m - 1), protein.atom(m - 1));
      }
    }
    }
  if ((!d_H.size()) || (!d_N.size())) {  // Checking only one should be enough
        throw Arguments::Exception("Selection is missing hydrogen or nitrogen atoms");
      }
  if ((!d_O.size()) || (!d_C.size())) {  // Checking only one should be enough
        throw Arguments::Exception("Selection is missing carbon or oxygen atoms");
      }

} //end Dssp::determineAtoms


void Dssp::calcHintra_init(utils::AtomSpecifier &protein)
{
  protein.sort();
  d_pbc = BoundaryParser::boundary(*d_sys, *d_args);  
  //this gather call does not do anything, 'cause we dont have coords...
  //d_pbc -> gather();
  for(unsigned int m=0; m<protein.size(); m++) {
    if(protein.name(m)=="CA") {
      d_CA.addAtom(protein.mol(m), protein.atom(m));
    }
  }
  if (!d_CA.size()) {  // Checking only one should be enough
        throw Arguments::Exception("Selection is missing alpha carbon atoms");
      }
}//end Dssp::calcHintra_init()

void Dssp::calcHb_Kabsch_Sander()
{
  acc_res.clear();
  don_res.clear();

  double rON=0, rCH=0, rOH=0, rCN=0;
  double q1=0.42, q2=0.20, f=33.2, E=0, cutoff=-0.5, distmin=0.05;

  Vec O(0.0,0.0,0.0);
  Vec C(0.0,0.0,0.0);
  for (int i = 0; i < (int) d_O.size(); ++i) {
    for (int j = 0; j < (int) d_H.size(); ++j) {
      if (d_O.mol(i) == d_H.mol(j) && d_O.mol(i) == d_N.mol(j) && d_O.mol(i) == d_C.mol(i)) {
    	O = d_pbc->nearestImage(*d_N.coord(j), *d_O.coord(i), d_sys->box());
    	rON = (*d_N.coord(j) - O).abs();
    	C = d_pbc->nearestImage(*d_H.coord(j), *d_C.coord(i), d_sys->box());
    	rCH = (*d_H.coord(j) - C).abs();
    	O = d_pbc->nearestImage(*d_H.coord(j), *d_O.coord(i), d_sys->box());
    	rOH = (*d_H.coord(j) - O).abs();
       	C = d_pbc->nearestImage(*d_N.coord(j), *d_C.coord(i), d_sys->box());
    	rCN = (*d_N.coord(j) - C).abs();
	    if (rON < distmin || rCH < distmin || rOH < distmin || rCN < distmin) {
	      E=-9.9;
	    }
	    else {
	      E = q1 * q2 * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN) * f;
        }

	if ((E < cutoff)
	    && (abs(
		d_O.resnum(i) + d_resOffSets[d_O.mol(i)] - d_H.resnum(j)
		    - d_resOffSets[d_H.mol(j)])) > 1) {
	  acc_res.push_back(d_O.resnum(i) + d_resOffSets[d_O.mol(i)]);
	  don_res.push_back(d_H.resnum(j) + d_resOffSets[d_H.mol(j)]);
	}
      }
    }
  }
} // end Dssp::calcHb_Kabsch_Sander()

void Dssp::calc_Helices()
{
  vector<int> turn3, turn4, turn5;
  vector<int> helix3_tmp, helix4_tmp, helix5_tmp;

  helix3.clear();
  helix4.clear();
  helix5.clear();
  Turn.clear();

  // loop over existing hydrogen bonds to identify elementary H-bond patterns
  // begin with three types of H-bonded turns
  // fill up Turn - may contain duplicates if one residue H-bonds to more than one 
  // other residue
  for (int i=0; i < (int) acc_res.size(); ++i) {
  
     if (don_res[i] == acc_res[i] + 3) {
       turn3.push_back(acc_res[i]);
       for (int j=1; j<3; j++) {
         Turn.push_back(acc_res[i]+j);
         }
     }
     if (don_res[i] == acc_res[i] + 4) {
       turn4.push_back(acc_res[i]);
       for (int j=1; j<4; j++) {
         Turn.push_back(acc_res[i]+j);
         }
     }
     if (don_res[i] == acc_res[i] + 5) {
       turn5.push_back(acc_res[i]);
       for (int j=1; j<5; j++) {
         Turn.push_back(acc_res[i]+j);
         }
     }
  }

  // see if turns form helices
  for (int i=0; i < (int) turn3.size()-1; ++i) {
    if (turn3[i] == (turn3[i+1] - 1)) {
      helix3_tmp.push_back(turn3[i]+1);
      helix3_tmp.push_back(turn3[i]+2);
      helix3_tmp.push_back(turn3[i]+3);
    }
  }
  for (int i=0; i < (int) turn4.size()-1; ++i) {
    if (turn4[i] == (turn4[i+1] - 1)) {
      helix4_tmp.push_back(turn4[i]+1);
      helix4_tmp.push_back(turn4[i]+2);
      helix4_tmp.push_back(turn4[i]+3);
      helix4_tmp.push_back(turn4[i]+4);
    }
  }

  for (int i=0; i < (int) turn5.size()-1; ++i) {
    if (turn5[i] == (turn5[i+1] - 1)) {
      helix5_tmp.push_back(turn5[i]+1);
      helix5_tmp.push_back(turn5[i]+2);
      helix5_tmp.push_back(turn5[i]+3);
      helix5_tmp.push_back(turn5[i]+4);
      helix5_tmp.push_back(turn5[i]+5);
    }
  }
  // remove "duplicates", this will also sort them
  for (int i = 0; i < numres; ++i) {
    if (helix3_tmp.size() > 0) {
      for (int j=0; j < (int) helix3_tmp.size(); ++j) {
	if (helix3_tmp[j] == d_resnum[i] ) {
	  helix3.push_back(helix3_tmp[j]);
	  break;
	}
      }
    }
    if (helix4_tmp.size() > 0) {
      for (int j=0; j < (int) helix4_tmp.size(); ++j) {
        if (helix4_tmp[j] == d_resnum[i]) {
	  helix4.push_back(helix4_tmp[j]);
	  break;
	}
      }
    }
    if (helix5_tmp.size() > 0) {
      for (int j=0; j < (int) helix5_tmp.size(); ++j) {
	if (helix5_tmp[j] == d_resnum[i] ) {
	  helix5.push_back(helix5_tmp[j]);
	  break;
	}
      }
    }
  }
} // end Dssp::calc_Helices()

void Dssp::calc_Betas()
{
  vector<int> p_bridge_tmp, ap_bridge_tmp, p_bridge_tmp2, ap_bridge_tmp2;
  vector<int> bridge_tmp, extended_tmp;

  bridge.clear();
  extended.clear();
  Beta.clear();

  // identify single beta bridges (parallel or antiparallel)
  // loop over Hbonds i and j (they are _not_ residue numbers)
  for (int i=0; i < (int) acc_res.size(); ++i) {
    for (int j=0; j < (int) acc_res.size(); ++j) {
      if ((don_res[i] == acc_res[j]) && (acc_res[i] == (don_res[j]-2))) {
	if (abs(acc_res[i]+1 - don_res[i]) > 2) {
	  p_bridge_tmp.push_back(acc_res[i]+1);
	  p_bridge_tmp.push_back(don_res[i]);
	}
      }
      if ((acc_res[i] == (don_res[j]-2)) && (don_res[i] == (acc_res[j]+2))) {
	if (abs(acc_res[i]+1 - don_res[i]-1) > 2) {
	  ap_bridge_tmp.push_back(acc_res[i]+1);
	  ap_bridge_tmp.push_back(don_res[i]-1);
	}
      }
      if ((don_res[i] == acc_res[j]) && (acc_res[i] == don_res[j])) {
	if (abs(don_res[i] - acc_res[i]) > 2) {
	  ap_bridge_tmp.push_back(don_res[i]);
	  ap_bridge_tmp.push_back(acc_res[i]);
	}
      }
    }
  }
  // remove "duplicates" also for the beta-bridges, this will also sort them
  for (int i=0; i < numres; ++i) {
    if (p_bridge_tmp.size() > 0 ) {  
      for (int j=0; j < (int) p_bridge_tmp.size(); ++j) {
	if (p_bridge_tmp[j] == d_resnum[i] ) {
	  p_bridge_tmp2.push_back(p_bridge_tmp[j]);
	  break;
	}
      }
    }
    if (ap_bridge_tmp.size() > 0 ) {
      for (int j=0; j < (int) ap_bridge_tmp.size(); ++j) {
	if (ap_bridge_tmp[j] == d_resnum[i] ) {
	  ap_bridge_tmp2.push_back(ap_bridge_tmp[j]);
	  break;
	}
      }
    }
  }
  // isolated bridge or extended strand?
  for (int i=0; i < (int) p_bridge_tmp2.size(); ++i) {    
    if (p_bridge_tmp2[i+1] == (p_bridge_tmp2[i] + 1)) {
      extended_tmp.push_back(p_bridge_tmp2[i]);
      extended_tmp.push_back(p_bridge_tmp2[i+1]);
    }
    else if ((p_bridge_tmp2[i] != (p_bridge_tmp2[i-1] +1)) || (i == 0)){
      bridge_tmp.push_back(p_bridge_tmp2[i]);
    }
  }
  for (int i=0; i < (int) ap_bridge_tmp2.size(); ++i) {
    if (ap_bridge_tmp2[i+1] == (ap_bridge_tmp2[i] + 1)) {
      extended_tmp.push_back(ap_bridge_tmp2[i]);
      extended_tmp.push_back(ap_bridge_tmp2[i+1]);
    }
    else if ((ap_bridge_tmp2[i] != (ap_bridge_tmp2[i-1] +1)) || (i == 0)) {
      bridge_tmp.push_back(ap_bridge_tmp2[i]);
    }
  }
  // remove duplicates, fill up Beta vector
  for (int i=0; i < numres; ++i) {
    if (extended_tmp.size() > 0 ) {
      for (int j=0; j < (int) extended_tmp.size(); ++j) {
	if (extended_tmp[j] == d_resnum[i] ) {
	  extended.push_back(extended_tmp[j]);
	  Beta.push_back(extended_tmp[j]);
	  break;
	}
      }
    }
    if (extended_tmp.size() > 0 ) {
      for (int j=0; j < (int) bridge_tmp.size(); ++j) {
	if (bridge_tmp[j] == d_resnum[i] ) {
	  bridge.push_back(bridge_tmp[j]);
	  Beta.push_back(bridge_tmp[j]);
	  break;
	}
      }
    }
  }
} //end Dssp::calc_Betas()

void Dssp::calc_Bends()
{
  gmath::Vec tmpA, tmpB; 
  double angle = 0.0;
  Bend.clear();
  d_pbc -> gather();

  for (int i = 2; i < (int) d_CA.size() - 2; ++i) {
    if (d_CA.mol(i - 2) == d_CA.mol(i) && d_CA.mol(i) == d_CA.mol(i + 2)) {
      tmpA = (*d_CA.coord(i - 2) - *d_CA.coord(i));
      tmpB = (*d_CA.coord(i + 2) - *d_CA.coord(i));
      angle = acos(
	  (tmpA.dot(tmpB)) / (sqrt(tmpA.dot(tmpA)) * (sqrt(tmpB.dot(tmpB)))))
	  * 180 / 3.1416;
      if (angle < 110) {
	Bend.push_back(d_resnum[i]);
      }
    }
  }
} //end Dssp::calc_Bends()

void Dssp::filter_SecStruct()
{
  vector<int>::iterator iter;

  turn.clear();

  // remove duplicates in Turn => turn
  for (int i=0; i < numres; ++i) {
    for (int j=0; j < (int) Turn.size(); ++j) {
      if (Turn[j] == d_resnum[i] ) {
	turn.push_back(d_resnum[i]);
	break;
      }
    }
  }
  
  // remove duplicates, priorities are h4>bridge>strand>h3>h5>turn>bend  
  // in the newest dssp version (2.1.0) by Maarten Hekkelman h5 (pi-helix) has been 
  // moved to the front -> not yet implemented here
  std::vector<std::vector<int>* > classes;
  classes.push_back(&helix4);
  classes.push_back(&bridge);
  classes.push_back(&extended);
  classes.push_back(&helix3);
  classes.push_back(&helix5);
  classes.push_back(&turn);
  classes.push_back(&Bend);
  
  for (unsigned int c = 0; c < classes.size()-1; c++) {
    for (unsigned int i = 0; i < classes[c]->size(); ++i) {
      for (unsigned int j = c+1; j < classes.size(); j++) {
        for (iter=classes[j]->begin(); iter!=classes[j]->end(); ++iter){
          if (*iter == (*classes[c])[i]) {
	        (*classes[j]).erase(iter);
	        --iter;
          }
        }
      }
    }
  }
} //end Dssp::filter_SecStruct()

void Dssp::keepStatistics()
{
  int typeIndex = 0;
  for (unsigned int i = 0; i < helix3.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == helix3[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;
  
  for (unsigned int i = 0; i < helix4.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == helix4[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;
  
  for (unsigned int i = 0; i < helix5.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == helix5[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;

  for(unsigned int i=0; i< turn.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == turn[i]) {
	index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;

  for (unsigned int i = 0; i < extended.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == extended[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;
  
  for (unsigned int i = 0; i < bridge.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == bridge[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;

  for (unsigned int i = 0; i < Bend.size(); ++i) {
    unsigned int index = 0;
    for (unsigned int z = 0; z < d_resnum.size(); z++) {
      if (d_resnum[z] == Bend[i]) {
        index = z;
      }
    }
    ++summary[index][typeIndex];
  }
  ++typeIndex;

  d_numFrames++;
  
}

void Dssp::writeToFiles(double time)
{
  for (int i=0; i < (int) helix3.size(); ++i) {
    timeseries3Helix << setw(10) << time << setw(10) << helix3[i]+1<< endl;
  }
  for (int i=0; i < (int) helix4.size(); ++i) {
    timeseries4Helix << setw(10) << time << setw(10) << helix4[i]+1<< endl;
  }
  for (int i=0; i < (int) helix5.size(); ++i) {
    timeseries5Helix << setw(10) << time << setw(10) << helix5[i]+1<< endl;
  }
  for (int i=0; i < (int) turn.size(); ++i) {
    timeseriesTurn << setw(10) << time << setw(10) << turn[i]+1<< endl;
  }
  for (int i=0; i < (int) extended.size(); ++i) {
    timeseriesBStrand << setw(10) << time << setw(10) << extended[i]+1 << endl;
  }
  for (int i=0; i < (int) bridge.size(); ++i) {
    timeseriesBBridge << setw(10) << time << setw(10) << bridge[i]+1 << endl;
  }
  for (int i=0; i < (int) Bend.size(); ++i) {
    timeseriesBend << setw(10) << time << setw(10) << Bend[i]+1 << endl;
  }
} //end Dssp::writeToFiles()

 void Dssp::opents(string fi1, string fi2, string fi3, string fi4, string fi5, string fi6, string fi7)
{
  timeseriesTurn.open(fi1.c_str());
  timeseries3Helix.open(fi2.c_str());
  timeseries4Helix.open(fi3.c_str());
  timeseries5Helix.open(fi4.c_str());
  timeseriesBBridge.open(fi5.c_str());
  timeseriesBStrand.open(fi6.c_str());
  timeseriesBend.open(fi7.c_str());
}

 void Dssp::readframe()
{
    InG96 icc;

    try{
      d_args -> check("ref",1);
      Arguments::const_iterator iterr=d_args -> lower_bound("ref");
      icc.open((iterr->second).c_str());
    }
    catch(const Arguments::Exception &){
      d_args -> check("traj",1);
      Arguments::const_iterator iterr=d_args -> lower_bound("traj");
      icc.open((iterr->second).c_str());
    }
    icc.select("ALL");
    icc >> *d_sys;
    icc.close();
    
}

Dssp::Dssp(gcore::System &sys, args::Arguments &args)
{
  d_sys=&sys;
  d_args=&args;
  d_O = AtomSpecifier(sys);
  d_H = AtomSpecifier(sys);
  d_C = AtomSpecifier(sys);
  d_N = AtomSpecifier(sys);  
  d_CA = AtomSpecifier(sys); 
  d_numFrames=0;
  
  opents("Turn.out", "3-Helix.out", "4-Helix.out", 
	 "5-Helix.out", "Beta-Bridge.out", "Beta-Strand.out", 
	 "Bend.out");
}

void Dssp::writeSummary(std::ostream & of)
{
  vector<int> average(7);
  of << "# Analysed " << d_numFrames << " structures\n#\n"
     << "#  res.      3-Helix      4-Helix      5-Helix         "
     << "Turn     B-Strand     B-Bridge         Bend\n"
     << "#            #     %      #     %      #     %      #  "
     << "   %      #     %      #     %      #     %\n";
  of.setf(ios::floatfield, ios::fixed);
  of.precision(1);
 
  for(unsigned int i=0; i< summary.size(); i++){
    of << setw(7) << d_resnum[i] + 1;
    for(unsigned int j=0; j < summary[i].size(); ++j){
      of << setw(7) << summary[i][j]
	 << setw(6) << 100*double(summary[i][j])/d_numFrames;
      average[j]+=summary[i][j];
    }
    of << endl;
  }
  of << "#\n";
  of << "#            3-Helix      4-Helix      5-Helix         "
     << "Turn     B-Strand     B-Bridge         Bend\n";
  of << "# protein     ";
  for(unsigned int i=0; i< average.size(); ++i)
    of << setw(6) << 100*double(average[i])/d_numFrames/summary.size()
       << "       ";
  of << endl;
}

void Dssp::calcnumres(utils::AtomSpecifier &protein, const System &sys)
{
  protein.sort();
  numres=0;
  d_resnum.clear();
  d_resOffSets.clear();
  d_resOffSets.push_back(0);
  int currentResNums = 0;

  for(int m = 0; m < sys.numMolecules() - 1; m++) {
    currentResNums += sys.mol(m).topology().resNum(sys.mol(m).numAtoms()-1) + 1;
    d_resOffSets.push_back(currentResNums);
  }

  for(unsigned int i=0, j=i+1; i<protein.size()-1; i++, j++) {

    while (protein.resnum(i) + d_resOffSets[protein.mol(i)]
	== protein.resnum(j) + d_resOffSets[protein.mol(j)]
	&& j < protein.size() - 1) {
      j++;
    }
    i=--j;
    numres++;
    d_resnum.push_back(protein.resnum(i) + d_resOffSets[protein.mol(i)]);
  }

  summary.resize(numres);
  for(unsigned int i=0; i< summary.size(); ++i){
    summary[i].resize(7);
  }
}


