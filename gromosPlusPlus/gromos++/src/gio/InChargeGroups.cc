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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InChargeGroups.cc
 * Author: bschroed
 * 
 * Created on March 8, 2018, 5:57 PM
 */
//#define DEBUG //FOR DEBUGING Uncomment
#include "InChargeGroups.h"

#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../gcore/AtomTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/System.h"
#include "../gromos/Exception.h"

using namespace std;
using namespace gcore;

namespace gio{
    //Constructor
    InChargeGroups::InChargeGroups() {
    }

    InChargeGroups::InChargeGroups(const InChargeGroups& orig){
    }

    InChargeGroups::InChargeGroups(string inChargeGroupFilePath) {
        this->inFilePath = inChargeGroupFilePath;   
        this->open(inChargeGroupFilePath);
        this->parseChargeGroupFile();
        
        //append all amino Acid residue names:
        //needed if terminus is directly an aa residue and has different atom set, than other aas.
        // APPEND resname of your terminal aa if not already here!         
        // TODO make an AA dicts
        
        aminoAcidSet.insert("ALA");
        aminoAcidSet.insert("LEU");
        aminoAcidSet.insert("GLY");
        aminoAcidSet.insert("SER");
        aminoAcidSet.insert("CYS");
        aminoAcidSet.insert("THR");
        aminoAcidSet.insert("ILE");
        aminoAcidSet.insert("VAL");
        aminoAcidSet.insert("VAL");
        aminoAcidSet.insert("ASN");
        aminoAcidSet.insert("GLN");
        aminoAcidSet.insert("PHE");
        aminoAcidSet.insert("LYS");
        aminoAcidSet.insert("TRP");
        aminoAcidSet.insert("LYN");
        aminoAcidSet.insert("ASH");
        aminoAcidSet.insert("ASP");
        aminoAcidSet.insert("MET");
        aminoAcidSet.insert("GLN");
        aminoAcidSet.insert("GLU");
        aminoAcidSet.insert("TYR");
        aminoAcidSet.insert("PRO");
        aminoAcidSet.insert("HIP");
        aminoAcidSet.insert("HID");
        aminoAcidSet.insert("HIE");
        aminoAcidSet.insert("CYM");
        aminoAcidSet.insert("CYX");        
        aminoAcidSet.insert("ARG");
    }

    InChargeGroups::~InChargeGroups() {

    }

    void InChargeGroups::parseChargeGroupFile(){
      if (!stream()) {
        ostringstream msg;
        msg << "Could not open AMBER topology file." << name();
        throw gromos::Exception("AmberTopology", msg.str());
      }

      // First read the whole file into the map
      vector<string> buffer;
      vector<string>::const_iterator it;
      
      while (!stream().eof()) {
        this->getblock(buffer);

        if (buffer.size()) {
          //std::cout << buffer[0] << std::endl; //comment
          d_blocks[buffer[0]] = buffer;
          buffer.clear();
        }
      }
      
      //translate blocks to columns
      map<string, map< int, pair<int, set<string>> > > resChargeGAtoms;
      for(auto& tmp : d_blocks){
            string resName = tmp.first;
            map< int, pair<int, set<string>>> chargeGroupAtoms;
            map< int, pair<int, set<string>> >::iterator it;
            #ifdef DEBUG
                cerr << "RES: " << resName<<endl; //comment
            #endif
            for(auto& row: tmp.second ){

                stringstream columns(row);
                string atomName;
                int chargeGroup;
                int chargeGroupCharge;
                
                //read Columns of block in charge groupfile
                columns >> atomName;
                columns >> chargeGroup;
                columns >> chargeGroupCharge;
                
                //ensure no title and no and END
                if(! atomName.compare("END") or ! atomName.compare(resName) ){             
                    continue;
                }
                
                //build up new structure: chargeGroupAtoms
                it = chargeGroupAtoms.find(chargeGroup);
                if(it == chargeGroupAtoms.end()){//charge group is not contained yet
                    set<string> atomNames;
                    atomNames.insert(atomName);
                    pair<int, set<string>> probs(chargeGroupCharge, atomNames);
                    chargeGroupAtoms.insert(pair<int, pair<int, set<string>>>(chargeGroup,probs));
                }
                else{
                  chargeGroupAtoms[chargeGroup].second.insert(atomName);
                }
            }
            resChargeGAtoms.insert(pair<string, map< int, pair<int, set<string>>>>(resName,chargeGroupAtoms));     
        }
        resChargeGroups = resChargeGAtoms;
    };
     
    map<int, pair<string,  double> >  InChargeGroups::spreadChargeInGroup(map<int, pair<string,  double> >  atomCharge, double chargeSum, int atomsNum, int resID, string resName, int chargeGroupID,  int ChargeGroupCharge /*= 0*/, set<string> exceptions /*= set<string>()*/){
        /*
         * spread charge diference of a charge group over all chargegroup atoms. 
         */
        
        double totalChargeDev = chargeSum - ChargeGroupCharge;
        map<int, pair<string,  double> > result;
        double chargeDiff = totalChargeDev / (atomsNum - exceptions.size()); 
        
        //testing
#ifdef DEBUG
        cerr <<"\tCorrect ChargeDev (chargeGroupCharge - chargeSum) : " <<totalChargeDev<<endl;//comment
        cerr <<"\tCorrect atoms with: " <<chargeDiff<<endl;//comment
#endif
        // Give warning if deviation in a group is very high
        double sumThreshold = 0.05;
        if (totalChargeDev > sumThreshold or totalChargeDev < 0-sumThreshold){
        	cerr << "WARNING:\t chargeDev is " << totalChargeDev << " (bigger than "<< sumThreshold<<" / "<< 0-sumThreshold<<") - for residue: "<< resID << " " << resName  <<"\t charge group: " << chargeGroupID << " \n\t\t The sum of all partial Charges should be: "<<ChargeGroupCharge << "\t BUT it is: "<< chargeSum <<endl<< endl;
        }
        
        for(map<int, pair<string,  double> >::iterator atomCit = atomCharge.begin(); atomCit != atomCharge.end(); atomCit++){
            pair<string, double> nameCharge = atomCit->second;
            if(exceptions.find(nameCharge.first) != exceptions.end()){ // do not change this charge!
                continue; //don"t change charge
            }
            else{ //spread charge
                nameCharge.second = nameCharge.second-chargeDiff;
            }
            result[atomCit->first] = nameCharge;
        }
        return result;
        
    };
    
    MoleculeTopology InChargeGroups::writeToMoleculeTopo (MoleculeTopology mTopo, map<int, pair<string,  double> > chargeGroupAtoms){
        // write new charges onto Mtopo
        for( map<int, pair<string,  double> >::iterator atomIt = chargeGroupAtoms.begin(); atomIt != chargeGroupAtoms.end();  ++atomIt ){
                           int tmpID= int(atomIt->first);
                            mTopo.atom(tmpID).setCharge(atomIt->second.second); 
                            if(++atomIt == chargeGroupAtoms.end()){
                                atomIt--;
                                mTopo.atom(tmpID).setChargeGroup(1);
                                #ifdef DEBUG
                                    cerr <<"\t"<< tmpID << "\t" << mTopo.atom(tmpID).name() <<"\t" << mTopo.atom(tmpID).charge() << "\t"<< mTopo.atom(tmpID).chargeGroup() << endl; //comment
                                #endif
                                break;
                            }    
                            else{
                                atomIt--;
                                mTopo.atom(tmpID).setChargeGroup(0);
                                #ifdef DEBUG
                                    cerr <<"\t" << tmpID << "\t" << mTopo.atom(tmpID).name() <<"\t" << mTopo.atom(tmpID).charge() << "\t"<< mTopo.atom(tmpID).chargeGroup() << endl; //comment
                                #endif
                            } 

                       }
        return mTopo;
    }
    
    void InChargeGroups::getChargeGroups(){};
    
    MoleculeTopology InChargeGroups::mapChargeGroupsOnAtoms(map< int, pair<int,set<string>>> AtomCharges, int resID, MoleculeTopology mTopo){
        int chargeGroupID = 0;
        string resName = mTopo.resName(resID);
        
        map< int, pair<int, set<string> > > chargeGroupAtomsScheme = AtomCharges;        // chargegroupID, (total chargeof group, atom names in group)
        set<string> atomNames = set<string>(pair<int, set<string>>(chargeGroupAtomsScheme.find(chargeGroupID)->second).second);
        int totalChargeOfChargeGroup = int(pair<int, set<string>>(chargeGroupAtomsScheme.find(chargeGroupID)->second).first);
        
        map<int, pair<string,  double> > chargeGroupAtoms; //Collect all atoms here
        
        //process vars
        double chargeSum = 0;
        int atomsNum=0;

        bool protein = false;
        //check If we have protein here?
        if(protein){
            //if yes add terminus:
            if(resID == 1){
            }
            else if(resID == mTopo.numRes()){
            }
        }
        
        //go through atom layer in res and collect atoms according to the Chargegroups
        //Assumption all atoms are in order! and present
        #ifdef DEBUG
            cerr << resName << "\n"; //comment
        #endif
        bool prot_first_cg;
        bool lastRes;
        
        if(resID == 0 and (aminoAcidSet.find(resName) != aminoAcidSet.end()) and chargeGroupID == 0 ){
            #ifdef DEBUG
                cerr << "\tFIRST Chargegroup of an aa- Protein start!" << endl;
            #endif
            atomNames.erase("H");
            atomNames.insert("H1");
            atomNames.insert("H2");
            atomNames.insert("H3");
            totalChargeOfChargeGroup = 1;
        }
        #ifdef DEBUG
            if(chargeGroupID<chargeGroupAtomsScheme.size()){
                cerr <<"\n"<< resName << "\tresID\t"<< resID << " of "<< mTopo.numRes()-1<< "\n"; //comment
                cerr <<  "\tcharge group: " << chargeGroupID << " of " << chargeGroupAtomsScheme.size()-1 << endl; //comment
                cerr << "\ttotal charge of group: "<< totalChargeOfChargeGroup << endl;
                cerr << "\tcontains:\t";
                for (string name : atomNames){
                    cerr << "  " << name ;
                }
                cerr << endl;
            }
            cerr<<"\n\tbefore chargecorrections:\n";
        #endif
        for(int a = 0; a < mTopo.numAtoms(); ++a){
           if(mTopo.resNum(a) == resID){//ensure correct residue
                  
               AtomTopology atom = mTopo.atom(a);
               string atmName = atom.name();
               int atomID = a;
               double atmCharge = atom.charge();
               #ifdef DEBUG
                cerr << "\t" << atomID << "\t" << atmName <<"\t" << atmCharge << "\t"<< atom.chargeGroup() << endl; //comment
               #endif
               atom.setChargeGroup(0);
               
               //find all atoms for charge group
               if(atomNames.find(atmName) != atomNames.end()){
                   chargeGroupAtoms[atomID] = pair<string, double>(atmName, atmCharge);
                   
                   chargeSum +=  atmCharge; // store atmCharge
                   atomsNum++; // add atom
                   
                   if(atomNames.size() <= atomsNum){
                       #ifdef DEBUG
                        cerr << "\n\tafter_correction" << endl; //comment
                       #endif
                       //alterCharge Groups
                       chargeGroupAtoms = spreadChargeInGroup(chargeGroupAtoms, chargeSum, atomsNum, resID, resName, chargeGroupID, totalChargeOfChargeGroup );
 
                       //write Back to Topology
                       atom.setChargeGroup(1);
                       mTopo = writeToMoleculeTopo(mTopo, chargeGroupAtoms);

		       if(chargeGroupID >= chargeGroupAtomsScheme.size()-1){
			break;
		       }
		       else{
			       //reset vars for next iterations
			       atomNames = set<string>(pair<int, set<string>>(chargeGroupAtomsScheme.find(++chargeGroupID)->second).second);
			       totalChargeOfChargeGroup = int(pair<int, set<string>>(chargeGroupAtomsScheme.find(chargeGroupID)->second).first);
			       
			       chargeGroupAtoms.clear();
			       chargeSum = 0.0;
			       atomsNum =0;
			       prot_first_cg = false;
			       
				//show current charge group
				#ifdef DEBUG
				    if(chargeGroupID<chargeGroupAtomsScheme.size()){
					cerr <<"\n"<< resName << "\tresID\t"<< resID << " of "<< mTopo.numRes()-1<< "\n"; //comment
					cerr <<  "\tcharge group: " << chargeGroupID << " of " << chargeGroupAtomsScheme.size()-1 << endl; //comment
					cerr << "\ttotal charge of group: "<< totalChargeOfChargeGroup << endl;
					cerr << "\tcontains:\t";
					for (string name : atomNames){
					    cerr << "  " << name ;
					}
					cerr << endl;
				    }
				#endif

			       //Identify if last chargegroup is part of an terminal Amino Acid? - if yes add terminal OXT
			       if(resID == mTopo.numRes()-1  and (aminoAcidSet.find(resName) != aminoAcidSet.end()) and chargeGroupAtomsScheme.size()-1 == chargeGroupID){
				   #ifdef DEBUG
				    cerr << "\tLAST CHARGEGroup of protein - Protterm!" << endl;
				   #endif 
				   atomNames.insert("OXT");
				   totalChargeOfChargeGroup = -1;
			       }
		      }
                   }                            
               }   
               else{
                   cerr << "MISSED ATOM "<< atomID << " " << atmName << " in residue: \t" << resID << " " << resName <<endl;//comment
                   ostringstream msg;
                   msg << "MISSED ATOM "<< atomID << " " << atmName << " in residue: \t" << resID << " " << resName << " !" <<endl;
                   throw gromos::Exception("ChargeGroupsIN-mapTo atoms", msg.str());
               }
           }
        }
        #ifdef DEBUG
		cerr << "\n \t\tEND CHARGE GROUPS\n " << endl;
                cerr.flush();
        #endif
	return mTopo;
    }
    
    gcore::System InChargeGroups::mapChargeGroupsOnResidues(System sys){
        int last_a =0; // temporary storage for efficency
        bool first = true;
        System resultSys;
        //traverse all atoms in topo
        for (int m = 0; m < sys.numMolecules(); ++m){
            MoleculeTopology mTopo = sys.mol(m).topology();
            //go through res layer
            for (int r = 0; r < mTopo.numRes() ; ++r){
                string resName = mTopo.resName(r);
                
                map<string, map< int, pair<int, set<string>> > >::iterator it; //iterator for Charge Groups
                it = resChargeGroups.find(resName);
              
                //Is there a new ChargeGroup for the residue?
                if(it == resChargeGroups.end()){
                    cerr << "WARNING:\tResidue " << resName <<" is not contained in ChargeGroup file" << endl; //comment
                    continue;
                }
                else{
                    #ifdef DEBUG
                    cerr << "Found Residue " << it->first << endl;
                    #endif
                    mTopo = mapChargeGroupsOnAtoms( it->second, r, mTopo);
                    #ifdef DEBUG
                    cerr << "mapped Charges " << it->first << endl;
                    #endif
                } 
            }
            Molecule mol(mTopo);
            resultSys.addMolecule(mol);
        }
        return resultSys;
    };    
}

// TODO ! make sure, that charge Group is SET! 
