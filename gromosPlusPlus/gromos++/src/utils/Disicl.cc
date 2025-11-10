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

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <utility>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>

#include "../gio/InTopology.h"
#include "../gio/Ginstream.h"
#include "../gcore/System.h"
#include "../utils/AtomSpecifier.h"
#include "../utils/PropertyContainer.h"
#include "../gromos/Exception.h"
#include "../utils/Disicl.h"

using utils::Dscl;


vector<string> &Dscl::split(const string &s, char delim, vector<string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
}


vector<string> Dscl::split(const string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


void Dscl::readLibrary(gio::Ginstream &lib) {  
  map<std::string, std::vector<std::string> > libContent;
  std::vector<std::string> buffer;
  
  // check for valid block structure and read lib file contents
  while(!lib.stream().eof()) {
    lib.getblock(buffer);
    if(!lib.stream().eof()) {
      libContent.insert(std::pair<string,std::vector<std::string> >(buffer[0],buffer));
    }    
  }
  readLibAng(libContent["DSCLANG"]);
  readLibReg(libContent["DSCLREG"]);
  readLibClass(libContent["DSCLCLASS"]);    
}


void Dscl::setPeriodic(int start, int end) {
    periodic.resize(2) ;
    periodic[0]=start;
    periodic[1]=end;
}

void Dscl::readLibAng(std::vector<std::string> &buffer) {
  int minshift=0;
  int maxshift=0;
  if (buffer.size()) { 
      string angName;
      std::vector<string> angAtoms(4);
      std::vector<int> angResShifts(4);   
      for(unsigned int i=1; i< buffer.size()-1; i++) {
	    std::istringstream linestream(buffer[i]);
	    linestream >> angName;
	    for (int iii = 0; iii < 4; ++iii) { 
	      linestream >> angAtoms[iii]; 
	      angAtomsUniq.push_back(angAtoms[iii]);
	    }
	    
	    for (int iii = 0; iii < 4; ++iii) { linestream >>  angResShifts[iii];}	  
	    
	           
	    libAngNames.push_back(angName);
	    libAngAtoms.push_back(angAtoms);
	    libAngResShifts.push_back(angResShifts);
	    
	    minshift = *std::min_element(angResShifts.begin(), angResShifts.end());
	    maxshift = *std::max_element(angResShifts.begin(), angResShifts.end());
	    if ( minshift < minResShift) 
	        minResShift = minshift;
	    if ( maxshift > maxResShift) 
	        maxResShift = maxshift;
	        
      } 
    std::sort (angAtomsUniq.begin(), angAtomsUniq.end());
    angAtomsUniq.erase( unique(angAtomsUniq.begin(), angAtomsUniq.end() ),
                        angAtomsUniq.end() ); 
    // cout << "min/maxresshift " << minResShift << " " << maxResShift << endl;
    numAng=libAngNames.size();
    
    // check for complex atom specifications
    for (unsigned int i=0; i< angAtomsUniq.size(); i++) {
      string &s=angAtomsUniq[i];
      multAtom ma;

      // this split function only works with single and 
      // not double quotes for the delimiter!
      vector<string> substrings=split(s,';');

      if (substrings.size() > 1) {
        ma.id=s;
        angAtomsUniq.erase(angAtomsUniq.begin()+i);
        // since I removed one entry I have to go back one step to point to the next entry
        i--; 
        ma.defaultAtom=substrings[0];
        
        // divide residue and atomnames
        for (unsigned int j=1; j<substrings.size(); j++) {
          vector<string> sstring=split(substrings[j],':');
          if (sstring.size()<2) {
            throw gromos::Exception("disicl","ERROR in library angle specifications!");
          }
          
          // there could be more than one residue for a given atomname
          vector<string> ssstring=split(sstring[0],',');
          if (ssstring.size() >1) {
            for (unsigned int k=0; k<ssstring.size(); k++) {
              ma.atoms.push_back(sstring[1]);
              ma.residues.push_back(ssstring[k]);
            }
          }
          else {
            ma.atoms.push_back(sstring[1]);
            ma.residues.push_back(sstring[0]);
          }
        }
        multAtoms.push_back(ma);
      }
      else {
        bool badformat=false;
        // the atom specifications should not contain ":" or ","
        // if they are simple atoms
        if (s.find(":")!=std::string::npos) {
          badformat=true;
        }
        if (s.find(",")!=std::string::npos) {
          badformat=true;
        }
        if (badformat) {
          throw gromos::Exception("disicl", 
	      "Bad atom name in library: "+ s);
	    }
      }
    }

    //some feedback to see if multi-atom input was read correctly
    cerr<< "# Atoms in angle definitions:\n#\t";
    for (unsigned int i=0; i< angAtomsUniq.size(); i++) {
      cerr << angAtomsUniq[i] << "  ";
    }
    cerr << "\n";
    if (multAtoms.size()>0) {
      //cerr << "# Complex atom specifications:\n";
      for (unsigned int j=0; j<multAtoms.size();j++) {
        cerr << "#\t" << multAtoms[j].id << " -> using ";
        for (unsigned int i=0; i< multAtoms[j].residues.size(); i++) {
          cerr << multAtoms[j].residues[i] <<":" << multAtoms[j].atoms[i] << ", ";
        }
        cerr << " for all other residues: "<< multAtoms[j].defaultAtom << endl;
      }
    }
  }
  else
      throw gromos::Exception("disicl", 
	     "No DSCLANG-block in library!");
} // end readLibAng

void Dscl::readLibReg(std::vector<std::string> &buffer) {
  if (buffer.size()) {   
      // cout << "# Reading DSCLREG\n";
      std::string regName;
      std::vector<double> regLimits(2*numAng); 
	  // loop over input lines
      for(unsigned int i=1; i< buffer.size()-1; i++) {
	    std::istringstream linestream(buffer[i]);
	    linestream >> regName;
	    bool outside=false;
	    for (unsigned int j = 0; j < 2*numAng; ++j) { 
	      linestream >> regLimits[j];
	      if (regLimits[j] < periodic[0] || regLimits[j] > periodic[1]) {
	          outside=true;
	      }
	    }
	    if (outside) {
  	      cerr << "WARNING: region "<< regName << "( ";
  	      for (unsigned int i=0; i < regLimits.size(); i++) {
  	        cerr << regLimits[i] << " ";
  	      }
  	      cerr  <<") outside periodic range ("
	         << periodic[0] << " to " << periodic[1] << ")\n";
	    }
	    libRegions.push_back(regLimits);
	    libRegionNames.push_back(regName);
	        
      }
      
    // check region limits
    for (unsigned int i=0; i< libRegions.size(); i++) {
      for (unsigned int j=0; j<numAng; j=j+2) {
        if (libRegions[i][j*2] >= libRegions[i][j*2+1]) { 
            cerr << "ERROR: something is wrong with the limits of region "<< 
                     libRegionNames[i]<<"\n";
            cerr << libRegions[i][j*2] << " >= " << libRegions[i][j*2+1] << "\n";
            throw gromos::Exception("disicl", "bad region" ); 
        }
      }
    }
    
      
    numReg=libRegions.size();
    
    // test for overlapping region limits
    for (unsigned int i=0; i < numReg; i++) {
      for (unsigned int j=i+1; j < numReg; j++) {
        unsigned int overlap=0;
        for (unsigned int ang=0; ang < numAng; ang++) {
          double min_i=libRegions[i][ang*2];
          double max_i=libRegions[i][ang*2+1];
          double min_j=libRegions[j][ang*2];
          double max_j=libRegions[j][ang*2+1];

          if ((min_i >= min_j && min_i < max_j) ||
              (max_i > min_j && max_i <= max_j) ||
              (min_i <= min_j && max_i >= max_j) )
            overlap++;
        }
        if (overlap == numAng) {
          stringstream errstring;
          errstring << libRegionNames[i] << " and " << libRegionNames[j] << "\n";
          for (unsigned int ang=0; ang < numAng; ang++) {
            errstring<<"|"<<libRegions[i][2*ang]<<" "<<libRegions[i][2*ang+1]<<"|";
          }
          errstring<<"\n";
          for (unsigned int ang=0; ang < numAng; ang++) {
            errstring<<"|"<<libRegions[j][2*ang]<<" "<<libRegions[j][2*ang+1]<<"|";
          }
          errstring<<"\n";
          cerr << "WARNING: Two regions are overlapping: " << errstring.str();
          cerr << "The first region in the list that matches will have precedence.\n";
          //throw gromos::Exception("disicl", "Two regions are overlapping: " + errstring.str());
          
        }
      }
    }
      
  }
  else
      throw gromos::Exception("disicl", 
	     "No DSCLREG-block in library!");
} // end readLibReg



void Dscl::readLibClass(std::vector<std::string> &buffer) {
  typedef multimap<string,vector<string> >::value_type mapTypeStr;    
  if(buffer.size()) {   
      string className, class1, class2, classShort, classDef;
      vector<string> classNames;
      
      int classcnt=1;
      for(unsigned int i=1; i< buffer.size()-1; i++) {
        std::istringstream linestream(buffer[i]);
	    linestream >> className >> class1 >> class2 >> classShort;
	    classDef=class1 + "-" + class2;
	    classNames.push_back(className);
	    classNames.push_back(classShort);

        // make a map of unique class names to be able to access 
        // name through shortname later
	    if (classNameMap.insert(map<string,string>::value_type(classShort,className)).second) {
	      classShortnUniq.push_back(classShort);
          classNameMap.insert(map<string,string>::value_type(classShort,className));
          classNumMap.insert(map<string,int>::value_type(classShort,classcnt));
          classcnt++;
        }
	    
	    if (libClassNames.insert(mapTypeStr(classDef,classNames)).second)
          libClassNames.insert(mapTypeStr(classDef,classNames));
        else
          throw gromos::Exception("disicl", "ERROR: Two classes have the same definition: " + classNames[0] + " and " + libClassNames[classDef][0] + ": " +classDef);
        classNames.clear();
      }

    classShortnUniq.push_back(unclassString);

    numClass=classShortnUniq.size();
  }
  else
      throw gromos::Exception("disicl", 
	     "No DSCLCLASS-block in library!");
} // end readLibClass

void Dscl::print_angles() {
  cout << "# ***ANGLES***" << endl;
  for (unsigned int i=0; i < numAng; i++) {
    cout << "# " << libAngNames[i] << " ";
    for (int j=0; j < 4; j++) { cout << libAngAtoms[i][j]  << " "; }
    for (int j=0; j < 4; j++) { cout << libAngResShifts[i][j] << " "; }
    cout << endl;
  }
}

void Dscl::print_regions() {
  cout << "# ***REGIONS***" << endl;
  for (unsigned int i=0; i < numReg ; i++) {
    cout  << "# " << libRegionNames[i];
    for(unsigned int j=0; j < numAng*2; j++) {
      cout << "\t" << libRegions[i][j];
    }
    cout << endl;
  }
}

void Dscl::print_classes() {
  cout << "# ***CLASSES***" << endl;
  for ( map<string, vector<string> >::const_iterator iter=libClassNames.begin();
        iter != libClassNames.end(); ++iter ) {
    cout  << "# " << iter->first << " -> " << (iter->second)[0] << "\t( " << (iter->second)[1] << " )"<< endl;
  }
}

void Dscl::writeHeader() {
  dihts.open(tsFile.c_str());
  dihts.setf(ios::floatfield, ios::fixed);
  dihts << "# ";
  for (int i = 0; i < 87; i++) { dihts << "-"; }
  dihts << "\n#" << setw(9) << "Time";
  dihts << setw(8) << "Residue";
  for(unsigned int i = 0; i < numAng; i++)
    dihts << setw(12) << libAngNames[i];
  dihts << "    # " << setw(12) << "Region" << setw(30) << "Class" << endl;
  dihts << "# ";
  for (int i = 0; i < 87; i++) { dihts << "-"; }
  dihts << "\n"; 
}


void Dscl::determineAtoms(utils::AtomSpecifier &atomSelection, gcore::System &sys) {
  atomSelection.sort();
  
  typedef map<string,vector<string> >::value_type mapTypeString;
  
  // store first residue number of every molecule
  // with count running through all molecules as in gromosXX
  molStartRes.push_back(0);
  int molStartCnt=0;
  
  for (int m=0; m < sys.numMolecules(); m++) {  
    molStartCnt+=sys.mol(m).topology().numRes();
    molStartRes.push_back(molStartCnt);
    map<string,vector<string> > atomExists;
    
    // exclude molecules with fewer than the necessary 
    // number of residues to calculate a region
    if (sys.mol(m).topology().numRes()>(maxResShift-minResShift)) {
      vector<string> atomList (sys.mol(m).topology().numRes(), "");
  
      for (unsigned int i=0; i<angAtomsUniq.size(); i++) {
        atomExists.insert(mapTypeString(angAtomsUniq[i],atomList));
      }
      for (unsigned int i=0; i<multAtoms.size(); i++) {
        atomExists.insert(mapTypeString(multAtoms[i].id,atomList));
      }
    } 
    vecAtomExists.push_back(atomExists);
    atomExists.clear();   
  }
  
  // store the maximal residue number for each molecule
  // so later I have to loop only until there
  maxResSel.resize(sys.numMolecules(), 0);
  int prevmol = -1;
  
  for (unsigned int i=0; i < atomSelection.size(); i++) {
    // do only if molecule has enough residues to calculate all angles
    if (sys.mol(atomSelection.mol(i)).topology().numRes()>(maxResShift-minResShift)) {
      bool found=false;
      // for single atom specs
      for (unsigned int j=0; j < angAtomsUniq.size(); j++) {
        if (atomSelection.name(i) == angAtomsUniq[j]) {
          vecAtomExists[atomSelection.mol(i)][angAtomsUniq[j]][atomSelection.resnum(i)]=atomSelection.name(i);
          found=true;
        }
      }
      // for multiatom specs (alternative atoms for different residues)
      for (unsigned int j=0; j < multAtoms.size(); j++) {
        multAtom &ma = multAtoms[j];
        bool resfound=false;
        for (unsigned int r=0; r < ma.residues.size(); r++) {
          if (atomSelection.resname(i) == ma.residues[r]) {
            resfound=true;
            if (atomSelection.name(i) == ma.atoms[r]) {
              vecAtomExists[atomSelection.mol(i)][ma.id][atomSelection.resnum(i)]=atomSelection.name(i);  
              found=true;
            }
            break;          
          }
        }
        if (!resfound) {
          // check for default
          if (atomSelection.name(i) == ma.defaultAtom) {
            vecAtomExists[atomSelection.mol(i)][ma.id][atomSelection.resnum(i)]=atomSelection.name(i);
            found=true;
          }        
        }
      }
      if (found) {
        if (atomSelection.mol(i) != prevmol) {
          prevmol=atomSelection.mol(i);
          molSel.push_back(prevmol);
        }
        if (maxResSel[atomSelection.mol(i)] < atomSelection.resnum(i)) {
          maxResSel[atomSelection.mol(i)] = atomSelection.resnum(i);
        }
      }
    }
    //else {
    //  cerr << "molecule " <<   atomSelection.mol(i)+1  << " has only " << 
    //       sys.mol(atomSelection.mol(i)).topology().numRes() << " residue(s) and will not be considered \n";
    //}
  }
  
  //find the highest atom number of each molecule
  maxAtom.resize(sys.numMolecules(), 0);
  for (vector<int>::iterator it=molSel.begin(), to=molSel.end();
       it != to; it++) {
    int m = *it;
    for (int i=0; i<sys.mol(m).topology().numAtoms(); i++) {
      if (sys.mol(m).topology().resNum(i) > maxResSel[m]) {
        maxAtom[m]=i;
        break;
      }
    }
  }
}   


void Dscl::getPropSpec(gcore::System &sys, bound::Boundary * pbc) {
  bool makeSpec;
  stringstream ss;
  
  fragments.push_back(new fragment);
  for (unsigned int i=0; i< numAng; i++) {
    fragments.back()->propStore.push_back(PropertyContainer (sys,pbc));
  }
  fragments.back()->propNames.resize(numAng);
  
  
  for (vector<int>::iterator it=molSel.begin(), to=molSel.end();
       it != to; it++) {
    int molnum = *it;
    map<string, vector<string> > &molecule = vecAtomExists[molnum];
    for (int i=0-minResShift; i <= maxResSel[molnum]-maxResShift; i++) {
      makeSpec=true;
      
      // for each dihedral read from library, test if all atoms of
      // all angles needed for classification exist
      for (unsigned int ang=0; ang < numAng; ang++) {
        for (unsigned int atom=0; atom < libAngAtoms[ang].size(); atom++) {
          string &atomname = libAngAtoms[ang][atom];
          int const &shiftedResNum = i+libAngResShifts[ang][atom];
          if ((molecule[atomname][shiftedResNum]).empty()) {
            makeSpec = false;
          }
        }
      }
      if (makeSpec) { 
        int resid = i+molStartRes[molnum];
        if (fragments.back()->resIds.size() != 0) {
          int lastmol = fragments.back()->molNums.back();
          int lastres = fragments.back()->resIds.back();
          int lastresi = fragments.back()->resNums.back();

          // start a new fragment if we start a new molecule or 
          // residues are not consecutive
          if (lastmol != molnum+1 || lastres != resid) {
            fragments.push_back(new fragment);
            for (unsigned int j=0; j< numAng; j++) {
              fragments.back()->propStore.push_back(PropertyContainer (sys,pbc));
            }
            fragments.back()->propNames.resize (numAng);
          }
        } 
        
        // put residue numbers for which propspecs are going to be created
        // in a vector, +1 because internal numbering starts at 0
        fragments.back()->resIds.push_back(resid+1);
        fragments.back()->resNums.push_back(i+1);
        fragments.back()->molNums.push_back(molnum+1);
        
        for (unsigned int ang=0; ang < numAng; ang++) {
          vector<int> shifts = libAngResShifts[ang];
          vector<string> atomnames = libAngAtoms[ang];
          // adding +1 because internal numbering of molecules and residues starts at 0
          ss << "t%" <<molnum+1<<":res(" << i+1+shifts[0] << ":" << molecule[atomnames[0]][i+shifts[0]] << ");" <<
          molnum+1<<":res(" << i+1+shifts[1] << ":" << molecule[atomnames[1]][i+shifts[1]] << ");" <<
          molnum+1<<":res(" << i+1+shifts[2] << ":" << molecule[atomnames[2]][i+shifts[2]] << ");" <<
          molnum+1<<":res(" << i+1+shifts[3] << ":" << molecule[atomnames[3]][i+shifts[3]] << ") " ;
          fragments.back()->propStore[ang].addSpecifier(ss.str());
          fragments.back()->propNames[ang].push_back(ss.str());
          ss.str ("");
        }
      }
    }
  }
  numFrags=fragments.size();
  numResTot=0;
  for (unsigned int f =0; f<numFrags; f++ ) {
    fragments[f]->numRes=fragments[f]->resIds.size();
    numResTot+=fragments[f]->numRes;
  }
  if (numResTot==0) {
      throw gromos::Exception("disicl", "No match for the given angle specifications.\n\t Maybe check your atom or residue names or the format of the angle specifications.");
  }
  getAtoms(sys);
}


void Dscl::getAtoms(gcore::System &sys) {
  for (unsigned int f=0; f<fragments.size(); f++) {
    for (unsigned int i=0; i<fragments[f]->numRes; i++) {
      AtomSpecifier as(sys);
      stringstream ss;
      ss << fragments[f]->molNums[i] << ":res(" << fragments[f]->resNums[i] << ":a)";
      as.addSpecifier(ss.str());
      fragments[f]->resAtoms.push_back(as);
      as.clear();
      ss.str ("");
    }
  }
}

void Dscl::calcDih() {
 for (unsigned int f=0; f<fragments.size(); f++) {
  for (unsigned int ang=0; ang < numAng; ang++) {
    fragments[f]->propStore[ang].calc();
  }  
 }
}


double Dscl::modulo(double value, int min, int max) {
  int diff;
  if (min<=max) 
    diff = max-min;
  else 
    throw gromos::Exception("disicl", "modulo min >= max!!");
  while (value > max)
    value -= diff;
  while (value <= min)
    value += diff;
  return value;
}


void Dscl::classifyRegions() {
 for (unsigned int f = 0; f<fragments.size(); f++) {
  // convert property-objects to Double to be able to compare them to 
  //library values and move to the given period  (default -180to180)
  fragments[f]->torsions.clear();
  fragments[f]->torsions.resize(numAng);
  for (unsigned int ang=0; ang < numAng; ang++) {
    for (unsigned int i=0; i < fragments[f]->numRes; i++) {
      fragments[f]->torsions[ang].push_back(modulo(fragments[f]->propStore[ang][i]->getValue().scalar(),periodic[0],periodic[1]));
    }
  }
  
  // compare dihedrals to library regions
  // and return the first region that fits  
  fragments[f]->regions.clear();
  vector<string> regionstmp(fragments[f]->numRes);
  string region;
  
  for (unsigned  int i=0; i < fragments[f]->numRes; i++) {
    regionstmp[i]=unclassString;
    for (unsigned int regcnt=0; regcnt < numReg; regcnt++) {
      vector<double> limits=libRegions[regcnt];  
      bool regcheck=true;
      for (unsigned int ang=0; ang < numAng; ang++) {
        if (!(fragments[f]->torsions[ang][i] >= limits[2*ang] && fragments[f]->torsions[ang][i] < limits[2*ang+1])) {
          regcheck=false;
          break;
        }
      }
      if (regcheck) {
        regionstmp[i]=libRegionNames[regcnt];
        break;
      }
    }
    fragments[f]->regions.push_back(regionstmp[i]);
  }
 }
}


void Dscl::classifyClasses(double const &time) {
 for (unsigned int f = 0; f<fragments.size(); f++) {
  fragments[f]->classes.clear();
  fragments[f]->bfactors.clear();
  string classname;
  // the last residue can never be classified because we need the +1 region
  for (unsigned int i=0; i < fragments[f]->numRes-1; i++) {
      if (libClassNames[fragments[f]->regions[i]+"-"+fragments[f]->regions[i+1]].size()) {
        classname=libClassNames[fragments[f]->regions[i]+"-"+fragments[f]->regions[i+1]][1];
      } else {
        classname=unclassString;
      }
      (*(classts[classname]))<< setw(10) << time << setw(10)<< fragments[f]->resIds[i] << "\n";
      fragments[f]->classes.push_back(classname);  
      fragments[f]->bfactors.push_back(classNumMap[classname]); 
  }
  // last residue
  fragments[f]->classes.push_back(unclassString); 
  fragments[f]->bfactors.push_back(0.0);
 }
}

void Dscl::getBfactorValues(gcore::System &sys) {
 for (unsigned int f = 0; f<fragments.size(); f++) {
  for (unsigned int i=0; i<fragments[f]->numRes; i++) {
    for (unsigned int j=0; j<fragments[f]->resAtoms[i].size(); j++) {
      sys.mol(fragments[f]->molNums[i]-1).setBfac(fragments[f]->resAtoms[i].atom(j),fragments[f]->bfactors[i]);
      
    }
  }
 }
}


void Dscl::writeDihTs(double const &time) {
 for (unsigned int f = 0; f<fragments.size(); f++) {
  for (unsigned int i=0; i< fragments[f]->numRes; i++) {
    dihts << setw(10) << std::setprecision(2) << time;
    dihts << setw(8) << fragments[f]->resIds[i];
    for (unsigned int ang=0; ang<numAng; ang++) {
      dihts << setw(12) << std::setprecision(4) << fragments[f]->torsions[ang][i];
    }
    dihts << "    # " << setw(12) << fragments[f]->regions[i] << setw(30) << fragments[f]->classes[i]<< endl;
  }
 }
}



void Dscl::initTimeseries() {
  string classname;
  for (unsigned int i=0; i<numClass; i++) { 
    classname = classShortnUniq[i];
    ofstream* of = new ofstream;
    classts.insert(map<string,ofstream* >::value_type(classname,of));
    classts[classname]->open(("class_"+classname+".dat").c_str());
    (*(classts[classname])) << "# " << (classname) << " - " << classNameMap[classname] << "\n";
    (*(classts[classname]))<<"#     time   residue\n";
  }
}


void Dscl::closeTimeseries() {
  string classname;
  for (unsigned int i=0; i<numClass; i++) {
    classname = classShortnUniq[i];
    classts[classname]->close();
    //deallocate
    delete classts[classname];
    classts[classname] = 0;
  }  
}

void Dscl::initSummary() {
  typedef multimap<string,vector<int> >::value_type mapTypeInt;
  vector<int> counts (numResTot, 0);
  for (unsigned int i=0; i < classShortnUniq.size(); i++) {
    summary.insert(mapTypeInt(classShortnUniq[i],counts));
  }
  summary.insert(mapTypeInt(unclassString,counts));
}

string Dscl::pdbTitle() {
  stringstream ss;
  ss << "DISICL out pdb with class information in the B-factor column\n";
  for (unsigned int i=0; i<classNumMap.size(); i++) {
    ss << classNumMap[classShortnUniq[i]] << " - " << classShortnUniq[i]  << ";  ";
    
    if ((i % 7) == 6) ss << endl;
  }
  return ss.str();
}

void Dscl::keepStatistics() {
int rescounter=0;
for (unsigned int f = 0; f<fragments.size(); f++) {
  for (unsigned int i=0; i<fragments[f]->numRes-1; i++) {
    ++summary[fragments[f]->classes[i]][rescounter];
    rescounter++;
  }
 }
}

void Dscl::writeStatistics(unsigned int  frameNum, bool do_tser) {
  vector<int> sum(classShortnUniq.size());
  int *entry;
  stats.open(statFile.c_str());
  stats.setf(ios::floatfield, ios::fixed);
  stats.precision(1);
  
  stats << "# Number of frames: " <<  frameNum << "\n";
  for (unsigned int j=0; j<classShortnUniq.size()-1; j++) {
    stats << "# "<< classShortnUniq[j] << " ... " << classNameMap[classShortnUniq[j]]<< "\n";
  }
  stats << "# " << unclassString<<" ... unclassified\n";
  stats << "# \n";
  
  stats <<"#" << setw(5) << "resi"<< setw(4) << "mol";
  for (unsigned int j=0; j<classShortnUniq.size(); j++) {
    stats << setw(13) << classShortnUniq[j];
  }
  stats << "\n";
  stats  << "#         ";
  for (unsigned int j=0; j<classShortnUniq.size(); j++) {
    stats << "      #     %";
  }
  stats <<"\n";
  
  // I exclude the last residue of each molecule, as I can not classify a 
  // class without the next residue's region anyways
  int rescounter=0;
  for (unsigned int f = 0; f<fragments.size(); f++) {
    for (unsigned int i=0; i<fragments[f]->numRes-1; i++) { 
      stats << setw(7) << fragments[f]->resIds[i] << setw(3) << fragments[f]->molNums[i];
      for (unsigned int j=0; j<classShortnUniq.size(); j++) {
        entry = &summary[classShortnUniq[j]][rescounter];
        stats << setw(7) << *entry << setw(6) << 100*double(*entry)/ frameNum;
        sum[j]+=*entry;
      }
      stats << endl;
      rescounter++;
    }
  }
    
  stats << "# \n";
  stats << "#         ";
  for (unsigned int j=0; j<classShortnUniq.size(); j++) {
    stats << setw(13) << classShortnUniq[j];
  }
  stats << "\n";
  stats << setw(7) << "# summary:";
  
  for (unsigned int j=0; j<sum.size(); j++) {
    stats << setw(7) << sum[j] << setw(6) << 100*double(sum[j])/ frameNum/(rescounter);
  }
  stats << "\n";
  
  // averages at the end of the dihedral timeseries file
  if (do_tser) {
    dihts << "# Averages over run: ";
    for (unsigned int i=0; i<numAng; i++) {
      dihts << libAngNames[i] << "(<average> <rmsd> <error estimate>)  " ;
    }
    dihts << "\n"; 
    for (unsigned int f=0; f<fragments.size(); f++) { 
    for (unsigned int i=0; i<fragments[f]->numRes; i++) {
      dihts << "# Residue " << setw(4) << fragments[f]->resIds[i] << ": ";
      for (unsigned int j=0; j<numAng; j++) {
        dihts << fragments[f]->propNames[j][i] << " ";
      }
      dihts << endl << "# averages:     ";
      for (unsigned int j=0; j<numAng; j++) {
        dihts << fragments[f]->propStore[j][i]->average() << " ";
      } 
      dihts << endl;
    }
    }
    dihts.close(); 
  }
}
  
void Dscl::writePdbColorLegend() {
  std::ofstream os;
  string pdbName="colorlegend.pdb";
  os.open(pdbName.c_str());
  os << "TITLE DISICL color legend\n";
  os << "REMARK residue names correspond to class names\n";
  os << "REMARK to be loaded into a visualization program together\n";
  os << "REMARK with the DISICL output pdbs and colored by b-factor\n";
  for (unsigned int i=0; i < classNumMap.size(); i++) {
    string resname=classShortnUniq[i].substr(0,3);
    string atomname;
    if (classShortnUniq[i].size() >3) atomname=classShortnUniq[i].substr(3,6);
    else atomname="X";
    os << "ATOM" << setw(7) << i+1 << setw(3) << atomname << setw(6) <<  resname
       << setw(6) << i+1 
       << setw(8) << i+1 << ".000   0.000   0.000  1.00" 
       << setw(3) << classNumMap[classShortnUniq[i]] << ".00\n";
  }
  os << "END\n";
  os.close();  
}
