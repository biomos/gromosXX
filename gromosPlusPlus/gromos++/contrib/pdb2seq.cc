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
 * @file pdb2seq.cc
 * Creates the building block sequence as well as the pdb2g96 library file
 * from a pdb file only
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pdb2seq
 * @section breas Creates the building block sequence as well as the pdb2g96 library file
 * @author @ref ae @ref bh
 * @date 18.11.2010
 *
 * Here goes the documentation...
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@pdb</td><td>&lt;pdb file&gt; </td></tr>
 * <tr><td>\@pH</dt><td>&lt;pH value of the simulation box&gt; </td></tr>
 * <tr><td>\@gff</td><td>&lt;gromos force field version&gt; </td></tr>
 * <tr><td>[\@select</dt><td>&lt;atoms to be read from PDB: \"ATOMS\" (standard), \"HETATM\' or \"ALL\"&gt;]
 * <tr><td>[\@head</dt><td>&ltbuilding block (sequence) of head group, e.g. NH3+&gt;] </td></tr>
 * * <tr><td>[\@tail</dt><td>&ltbuilding block (sequence) of tail group, e.g. COO-&gt;] </td></tr>
 * </table>
 *
 *
 * Example using a specification file:
 * @verbatim
  cry
    @pdb     protein.pdb
    @gff     53a6
    @pH      7
 @endverbatim
 *
 * <hr>
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdio>

#include "../src/args/Arguments.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/AminoAcid.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InPDB.h"


// REMEMBER
// ========
//
// - check the constants below
// - adapt the InPDB to read also gromos PDBs


// CONSTANTS
// =========
//
const double SS_cutoff = 2.2;
/*
PDB avarage distance for the S-S bond is 2.03 Ansgtrom
Plus arbitrary 0.17 to account for deviations
Engh. et al. Int. Tables for Cryst. (2006), F, 18.3, pp382
*/
const double Hbond_dist = 3.5; //<<<< Check!!!


 //FUNCTION DECLARATIONS
 // =====================
 //
std::vector<std::string> findSS(gio::InPDB &myPDB);
std::vector<std::string> AcidOrBase(std::vector<std::string> seq, double pH,
        utils::gromosAminoAcidLibrary &gaal);
std::vector<std::string> EndGroups(gio::InPDB &myPDB, std::vector<std::string> seq, double pH,
        utils::gromosAminoAcidLibrary &gaal, std::string head, std::string tail);
std::vector<std::string> Histidine(gio::InPDB &myPDB, std::vector<std::string> seq, double pH,
        utils::gromosAminoAcidLibrary &gaal);
void writeResSeq(std::ostream &os, std::vector<std::string> seq);
void writePDB(gio::InPDB &myPDB, std::vector<std::string> seq);
void checkAdaptPDB(gio::InPDB &myPDB, utils::gromosAminoAcidLibrary &gaal);


using namespace std;
using namespace args;
using namespace gio;
using namespace utils;


int main(int argc, char **argv) {


  // DEFINE THE COMMAND LINE ARGUMENTS AND USAGE STRING
  // ==================================================
  //
  Argument_List knowns;
  knowns << "pdb" << "gff" << "pH" << "aalib" << "select" << "head" << "tail";
  //
  string usage = "# " + string(argv[0]);
  usage += "\n\t@pdb      <pdb file to be read>\n";
  usage += "\t@pH       <specification file for the symmetry transformations]>\n";
  usage += "\t[@select  <atoms to be read from PDB: \"ATOMS\" (standard), \"HETATM\' or \"ALL\">]\n";
  usage += "\t[@aalib   <amino acid library file>]\n";
  usage += "\t[@gff     <GROMOS force field version (if no @aalib sepcified): 45A4, 53A6>]\n";
  usage += "\t[@head    [<building block (sequence) of head group, e.g. NH3+>]]\n";
  usage += "\t[@tail    [<building block (sequence) of tail group, e.g. COO->]]\n";

  try {

    // READ THE COMMAND LINE ARGUMENTS
    // ===============================
    //
    Arguments args(argc, argv, knowns, usage);
	// this program is under development, so test if the @develop flag is there
	args.underDevelopment();
    //
    // which PDB file to be used (file name)?
    if (args.count("pdb") != 1) {
      throw gromos::Exception(argv[0], "specify exactly one pdb file (@pdb)");
    }
    //
    // the agromos force fiel version to be used
    string gffversion;
    if (args.count("gff") == 1) {
      gffversion = args.find("gff")->second;
    } else {
      throw gromos::Exception(argv[0], "specify exactly one GROMOS force field version (@gff)");
    }
    //
    // simulation intended to run at which pH value?
    double pH;
    {
      stringstream ss, css;
      if (args.count("pH") == 1) {
        ss.str(args.find("pH")->second);
        css.str(ss.str());
        ss >> pH;
        if (ss.fail() || ss.bad()) {
          stringstream msg;
          msg << "could not convert " << css.str() << " to a valid pH value";
          throw gromos::Exception(argv[0], msg.str());
        }
      } else {
        throw gromos::Exception(argv[0], "no or more than one value indicated as pH (@pH)");
      }
    }

    // selection of atoms read from PDB
    string select = "ATOM";
    if(args.count("select") > 0) {
      select = args["select"];
      if(select != "ATOM" && select != "HETATM" && select != "ALL") {
        stringstream msg;
        msg << select << " is not a proper selection of atoms to be read from pdb"
                " (@select), allowed is \"ATTOM\", \"HETATM\" or \"ALL\"";
        throw gromos::Exception(argv[0], msg.str());
      }
    }
    
    //HEAD or TAIL definitions
    string head = "noHead";
    string tail = "noTail";

    if(args.count("head") == 0) {
      head = "NHX";
    } else if(args.count("head") == 1) {
      head = args.find("head")->second;
    }else if(args.count("head")>1) {
      stringstream ss;
      for(Arguments::const_iterator it = args.lower_bound("head");
              it != args.upper_bound("head"); ++it) {
        ss << it->second << " ";
      }
      head = ss.str().substr(0,ss.str().size()-1);
    }

    if(args.count("tail") == 0) {
      tail = "COOX";
      //cerr << "new tail = " << tail << endl;
    } else if(args.count("tail") == 1) {
      tail = args.find("tail")->second;
    }else if(args.count("tail")>1) {
      stringstream ss;
      for(Arguments::const_iterator it = args.lower_bound("tail");
              it != args.upper_bound("tail"); ++it) {
        ss << it->second << " ";
      }
      tail = ss.str().substr(0,ss.str().size()-1);
    }

    // READ THE PDB FILE
    // =================
    //
    InPDB ipdb(args["pdb"]);
    ipdb.select(select);
    ipdb.read();
    ipdb.renumberRes();

    // BUILD/READ THE LIBRARY FILES
    // ============================
    //
    utils::gromosAminoAcidLibrary gaal;
    gaal.loadHardcoded45A4();
    //gaal.writeLibrary(cout, "  45A4 AMINO ACID LIBRARY");

    // RESIDUE SEQUENCE FROM PDB
    // =========================
    //
    // check if all residue name of the PDB are recognized in the AminoAcid library
    checkAdaptPDB(ipdb, gaal);
    // extract the PDB residue sequence
    vector<string> resSeq = ipdb.getResSeq();
    // check for disulfide briges
    resSeq = findSS(ipdb);
    // adapt the protonation state of the residues
    resSeq = AcidOrBase(resSeq, pH, gaal);
    
    // HISTIDINS
    // =============
    //
    // decide about the His protonation state
    // If HIS is base, then HISX can be HISA or HISB

    resSeq = Histidine(ipdb, resSeq, pH, gaal);

    writePDB(ipdb, resSeq);

    // add the head and tail groups (if specified in the input file to do so)
    //cerr << "head = " << head << endl << "tail = " << tail << endl;
    resSeq = EndGroups(ipdb, resSeq, pH, gaal, head, tail);
    
    // PRINT OUT ALL WARNINGS/ERRORS
    // =============================
    //
    //

    // WRITE THE SEQUENCE AND LIBRARIES
    // ================================
    //
    // - residue sequence
    // - pdb2g96 library
    // - "corrected" PDB
    //
    writeResSeq(cout, resSeq);

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



//FUNCTION DEFINITIONS
/**
 * Check for S-S bridges
 */
std::vector<std::string> findSS(InPDB &myPDB) {

  int num = 0;

  stringstream ss;

  vector<string> sequence = myPDB.getResSeq();
  for (unsigned int i = 0; i < myPDB.numAtoms(); ++i) {
    if (sequence[myPDB.getResNumber(i) - 1] == "CYS") {
      if (myPDB.getAtomName(i) == "SG") {
        gmath::Vec first_cys = myPDB.getAtomPos(i);
        for (unsigned int j = i + 1; j < myPDB.numAtoms(); ++j) {
          if (sequence[myPDB.getResNumber(j) - 1] == "CYS") {
            if (myPDB.getAtomName(j) == "SG") {
              gmath::Vec second_cys = myPDB.getAtomPos(j);
              double dist = (first_cys - second_cys).abs();
              if (dist < SS_cutoff) {
                sequence[myPDB.getResNumber(i) - 1] = "CYS1";
                sequence[myPDB.getResNumber(j) - 1] = "CYS2";
                ss << myPDB.getResNumber(i) << "-" << myPDB.getResNumber(j) << " ";
                num++;
              }
            }
          }
        }
      }
    }
  }
  
  if (ss.str() != ""){
    ofstream fout("ss.dat");
    fout << ss.str();
    fout.close();
  }


  //cout << num << " SS-briges found\n";

  return sequence;
}

vector<string> AcidOrBase(vector<string> seq, double pH, utils::gromosAminoAcidLibrary &gaal) {
  for(unsigned int i = 0; i < seq.size(); ++i) {
    if(seq[i] != "CYS1" && seq[i] != "CYS2") {
      double pKc = gaal.pKc(seq[i]);
      if(pKc > 0 && pH > pKc) {
        seq[i] = gaal.pdb2base(seq[i]);
      } else {
        seq[i] = gaal.pdb2acid(seq[i]);
      }
    }
  }
  return seq;
}

vector<string> EndGroups(InPDB &myPDB, vector<string> seq, double pH, 
        utils::gromosAminoAcidLibrary &gaal, string head, string tail){
  

  //Finding where to put end groups
  // By default, before the first residue and after the last residue
  vector<unsigned int> startposition;
  vector<unsigned int> endposition;
  if (head != "noHead") {
    startposition.push_back(1);
  }
  for(unsigned int i = 0; i < myPDB.numAtoms()-1; ++i) {
    if (myPDB.getChain(i) != myPDB.getChain(i + 1)) {
      if (head != "noHead") {
        endposition.push_back(myPDB.getResNumber(i));
      }
      if (tail != "noTail") {
        startposition.push_back(myPDB.getResNumber(i + 1));
      }
    }
  }
  // add the default end group
  if (tail != "noTail") {
    endposition.push_back(seq.size());
  }
  
  vector<string> start;
  vector<string> end;

  for (unsigned int i = 0; i<startposition.size(); ++i){
    if(head == "NHX"){
      double pKb = gaal.pKb(myPDB.getResName(startposition[i]-1));
      if(pH > pKb){
        start.push_back("NH2");
      }else{
        start.push_back("NH3+");
      }
    }else{
      start.push_back(head);
    }
  }
  for (unsigned int i = 0; i<endposition.size(); ++i){
    if(tail == "COOX"){
      double pKa = gaal.pKa(myPDB.getResName(endposition[i]-1));
      if(pH > pKa){
        end.push_back("COO-");
      }else{
        end.push_back("COOH");
      }
    }else{
      end.push_back(tail);
    }
  }


  vector<string> newSeq;
  unsigned int j = 0;
  unsigned int k = 0;

  
  //bool dirtyfix = false;
  for (unsigned int i = 0; i < seq.size(); i++) {
    //cout << "k: "<< k << " i: " << i << " j: "<< j << " endpos: "
    //          <<endposition[k]<<"startpos: " << startposition[j]<< endl;

    /*
    if (i == startposition[j] - 1) {
      newSeq.push_back(start[j]);
      newSeq.push_back(seq[i]);
      if (j < startposition.size()) {
        j++;
      }
    }
    
     

    if (i == endposition[k]-1) {
      newSeq.push_back(seq[i]);
      newSeq.push_back(end[k]);
      if (k < endposition.size()-1) {
        k++;
        dirtyfix = true;
      }
    }
    

   // cout << "k: "<< k << " i: " << i << " j: "<< j <<endl;

    //if ((i != endposition[k]-1  ) && (i != startposition[j-1] )) {
    if ((dirtyfix == false  ) && (i != startposition[j-1]-1 )) {
      newSeq.push_back(seq[i]);
      //cout << "k: "<< k << " i: " << i << " j: "<< j << " endpos: "
       //       <<endposition[k]<<"startpos: " << startposition[j]<< endl;
    }
    if ((dirtyfix == true  ) && (i != endposition[k-1]  )&& (i != startposition[j-1] )) {
      newSeq.push_back(seq[i]);
    }*/

    if (startposition.size() > 0 && i == startposition[j] - 1) {
      newSeq.push_back(start[j]);
      newSeq.push_back(seq[i]);
      if (j < startposition.size()) {
        j++;
      }
    } else if (endposition.size() > 0 && i == endposition[k]-1) {
      newSeq.push_back(seq[i]);
      newSeq.push_back(end[k]);
      if (k < endposition.size()-1) {
        k++;
        //dirtyfix = true;
      }
    } else {
      newSeq.push_back(seq[i]);
    }
  }
  return newSeq;

}

std::vector<std::string> Histidine(gio::InPDB &myPDB, std::vector<std::string> seq, double pH,
        utils::gromosAminoAcidLibrary &gaal){

  // Information is stored in this map
  //std::map<int, int> histypes;
  //histypes.clear();
  // The first integer is the residue number
  // The second integer is the code that specifies the state:
  /*
   * 0 - Nothing is close by (initial state)
   * 1 - ND1 has donor close by
   * 2 - ND1 has acceptor close by
   * 3 - NE2 has donor close by
   * 4 - NE2 has acceptor close by
   * 5 - both ND1 and NE2 have donor close by
   * 6 - both ND1 and NE2 have acceptor close by
   * 7 - ND1 has donor close by and NE2 has acceptor close by
   * 8 - NE2 has donor close by ND1 has acceptor close by
   */

  bool ND1_hasD = false;
  bool ND1_hasA = false;
  bool NE2_hasD = false;
  bool NE2_hasA = false;

  double donor_dist_1 = 9999;
  double acceptor_dist_1 = 9999;
  double donor_dist_2 = 9999;
  double acceptor_dist_2 = 9999;
  //histypes.insert(pair<int, int> (myPDB.getResNumber(i), 0));

  //debug
  //bool found_something = false;
  // found what?
  //vector<string> foundwhat;
  //vector<string> foundinres;
  //vector<double> foundwithdist;

  for(unsigned int i = 0; i < myPDB.numAtoms(); ++i) {
    if(myPDB.getResName(i) == "HIS") {
      //cout << "have I been here HISX" << endl;
      //if(histypes.empty()){
      //  histypes.insert(pair<int, int> (myPDB.getResNumber(i), 0));
      //}
      //if(myPDB.getAtomName(i) == "N" ){
      //  cout << myPDB.getResNumber(i) << "HIS"<<endl;
      //}

      if(myPDB.getAtomName(i) == "ND1") {
        //cout << "have I been here?" << endl;
        for(unsigned int j = 0; j < myPDB.numAtoms(); ++j) {
          if(myPDB.getResNumber(i) != myPDB.getResNumber(j)) {
            for(unsigned int k = 0; k < gaal.rHdonors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1]).size(); ++k) {
              //cout << "level ND1 after donor k-loop" << endl;
              if(myPDB.getAtomName(j) == gaal.rHdonors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1])[k]){
                double dist;
                dist = (myPDB.getAtomPos(i)-myPDB.getAtomPos(j)).abs();
                //found_something =true;
                //foundwhat.push_back(myPDB.getAtomName(j));
                //foundinres.push_back(myPDB.getResName(j));
                //foundwithdist.push_back(dist);
                if (dist < Hbond_dist && dist < donor_dist_1){
                  //cout << myPDB.getResNumber(i) <<"  ND1 check donor :" << dist
                  //        << " to " << myPDB.getResName(j) << myPDB.getResNumber(j) << endl;
                  donor_dist_1 = dist;
                  ND1_hasD = true;
                  //histypes.clear();
                  //histypes.insert(pair<int, int> (myPDB.getResNumber(i), 1));
                }
              }
            }
          }
        }
        for(unsigned int j = 0; j < myPDB.numAtoms(); ++j) {
          if(myPDB.getResNumber(i) != myPDB.getResNumber(j)) {
            for(unsigned int k = 0; k < gaal.rHacceptors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1]).size(); ++k) {
              //cout << "level ND1 after acceptor k-loop" << endl;
              if(myPDB.getAtomName(j) == gaal.rHacceptors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1])[k]){
                double dist;
                dist = (myPDB.getAtomPos(i)-myPDB.getAtomPos(j)).abs();
                //ound_something =true;
                //foundwhat.push_back(myPDB.getAtomName(j));
                //foundinres.push_back(myPDB.getResName(j));
                //foundwithdist.push_back(dist);
                if (dist < Hbond_dist && dist < donor_dist_1 && dist < acceptor_dist_1){
                  //cout << myPDB.getResNumber(i) <<"  ND1 check acceptor :" << dist
                  //        << " to " << myPDB.getResName(j) << myPDB.getResNumber(j) << endl;
                  acceptor_dist_1 = dist;
                  ND1_hasA = true;
                  ND1_hasD = false;
                  //histypes.clear();
                  //histypes.insert(pair<int, int> (myPDB.getResNumber(i), 2));
                }
              }
            }
          }
        }       
      }
      if(myPDB.getAtomName(i) == "NE2") {
        for(unsigned int j = 0; j < myPDB.numAtoms(); ++j) {
          if(myPDB.getResNumber(i) != myPDB.getResNumber(j)) {
            for(unsigned int k = 0; k < gaal.rHdonors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1]).size(); ++k) {
              //cout << "level NE2 after donor k-loop" << endl;
              if(myPDB.getAtomName(j) == gaal.rHdonors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1])[k]){
                double dist;
                dist = (myPDB.getAtomPos(i)-myPDB.getAtomPos(j)).abs();
                //found_something =true;
                //foundwhat.push_back(myPDB.getAtomName(j));
                //foundinres.push_back(myPDB.getResName(j));
                //foundwithdist.push_back(dist);
                if (dist < Hbond_dist && dist < donor_dist_2){
                  NE2_hasD = true;
                  //cout << myPDB.getResNumber(i) <<"  NE2 check donor :" << dist
                  //        << " to " << myPDB.getResName(j) << myPDB.getResNumber(j) << endl;
                  donor_dist_2 = dist;
                  //int a = 0;
                  //if(histypes.find(myPDB.getResNumber(i))->second == 1){
                  // a = 5;
                  //} else if (histypes.find(myPDB.getResNumber(i))->second == 2){
                  //  a = 8;
                  //} else if (histypes.find(myPDB.getResNumber(i))->second != 1 && histypes.find(myPDB.getResNumber(i))->second != 2){
                  //  a = 3;
                  //}
                  //histypes.clear();
                  //histypes.insert(pair<int, int> (myPDB.getResNumber(i), a));
                }
              }
            }
          }
        }
        for(unsigned int j = 0; j < myPDB.numAtoms(); ++j) {
          if(myPDB.getResNumber(i) != myPDB.getResNumber(j)) {
            for(unsigned int k = 0; k < gaal.rHacceptors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1]).size(); ++k) {
              //cout << "level NE2 after acceptor k-loop" << endl;
              if(myPDB.getAtomName(j) == gaal.rHacceptors(myPDB.getResName(j), seq[myPDB.getResNumber(j)-1])[k]){
                double dist;
                dist = (myPDB.getAtomPos(i)-myPDB.getAtomPos(j)).abs();
                //found_something =true;
                //foundwhat.push_back(myPDB.getAtomName(j));
                //foundinres.push_back(myPDB.getResName(j));
                //foundwithdist.push_back(dist);
                if (dist < Hbond_dist){
                  //cout << myPDB.getResNumber(i) << " with " << myPDB.getResName(j)<< myPDB.getResNumber(j)
                  //        << " with dist: " << dist << " and DD2 "<< donor_dist_2 << " AD2 " << acceptor_dist_2 << endl;
                }
                if (dist < Hbond_dist && dist < donor_dist_2 && dist < acceptor_dist_2){
                  //cout << myPDB.getResNumber(i) <<"  NE2 check acceptor :" << dist
                  //        << " to " << myPDB.getResName(j) << myPDB.getResNumber(j) << endl;
                  acceptor_dist_2 = dist;
                  NE2_hasA = true;
                  NE2_hasD = false;
                  //double a;
                  //if(histypes.find(myPDB.getResNumber(i))->second == 1){
                  //  a = 7;
                  //} else if (histypes.find(myPDB.getResNumber(i))->second == 2){
                  //  a = 6;
                  //} else if (histypes.find(myPDB.getResNumber(i))->second != 1 && histypes.find(myPDB.getResNumber(i))->second != 2){
                  //  a = 4;
                  //}
                  //histypes.clear();
                  //histypes.insert(pair<int, int> (myPDB.getResNumber(i), a));
                }
              }
            }
          }
        }       
      }
      
      //int histype = histypes.find(myPDB.getResNumber(i))->second;
      //histypes.clear();
      //cout << histype << " << histype" << endl;

      if(i == myPDB.numAtoms() - 1 || myPDB.getResNumber(i) != myPDB.getResNumber(i + 1)) {
        /*switch(histype) {
          case 0:
          {
            //BY DEFAULT: put HISB as we have no clue! :)
            seq[myPDB.getResNumber(i) - 1] = "HISB";
            break;
          }
          case 1:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISB";//Fine
            break;
          }
          case 2:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISA";//Fine
            break;
          }
          case 3:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISA";//Fine
            break;
          }
          case 4:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISB";//Fine
            break;
          }
          case 5:
          {
            //BY DEFAULT: put HISB as we have no clue! :)
            if(donor_dist_1 < donor_dist_2) {
              seq[myPDB.getResNumber(i) - 1] = "HISB";
            } else {
              seq[myPDB.getResNumber(i) - 1] = "HISA";
            }
            break;
          }
          case 6:
          {
            //BY DEFAULT: put HISB as we have no clue! :)
            seq[myPDB.getResNumber(i) - 1] = "HISB";
            break;
          }
          case 7:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISB";
            break;
          }
          case 8:
          {
            seq[myPDB.getResNumber(i) - 1] = "HISA";
            break;
          }
          default:
          {
            break;
          }
        }*/

        /*
         * Testing the conditions
         * - ND1 can only have D or A (see condition above)
         * - NE2 can only have D or A (see condition above)
         *
         */

        if (ND1_hasD){
          if (NE2_hasD){
            if (donor_dist_1 <= donor_dist_2){
              seq[myPDB.getResNumber(i) - 1] = "HISB";
            }else{
              seq[myPDB.getResNumber(i) - 1] = "HISA";
            }
          }else if (NE2_hasA){
            seq[myPDB.getResNumber(i) - 1] = "HISB"; //definitely!!!
          }else{
            seq[myPDB.getResNumber(i) - 1] = "HISB";
          }
        }else if (ND1_hasA){
          if(NE2_hasD) {
            seq[myPDB.getResNumber(i) - 1] = "HISA"; //definitely!!!
          } else if(NE2_hasA) {
            if(acceptor_dist_1 <= acceptor_dist_2) {
              seq[myPDB.getResNumber(i) - 1] = "HISA";
            } else {
              seq[myPDB.getResNumber(i) - 1] = "HISB";
            }
          } else {
            seq[myPDB.getResNumber(i) - 1] = "HISA";
          }
        } else {
          if(NE2_hasD) {
            seq[myPDB.getResNumber(i) - 1] = "HISA";
          } else if(NE2_hasA) {
            seq[myPDB.getResNumber(i) - 1] = "HISB";
          } else {
            seq[myPDB.getResNumber(i) - 1] = "HISB"; // by definition
          }
        }

        //Resetting the distances:
        donor_dist_1 = 9999;
        acceptor_dist_1 = 9999;
        donor_dist_2 = 9999;
        acceptor_dist_2 = 9999;

        //cout << ND1_hasD << " " << ND1_hasA << " " << NE2_hasD << " " << NE2_hasA
        //        << myPDB.getResName(i) << myPDB.getResNumber(i) << " found: " << found_something << endl;

        //if(found_something) {
        //  for(unsigned int i = 0; i < foundwhat.size(); i++) {
        //    if(foundwithdist[i] < 3.5) {
        //      cout << foundwhat[i] << "-" << foundinres[i] << "-> " << foundwithdist[i] << endl;
        //    }
        //  }
        //  cout << "Found this number of stuff   :" << foundwhat.size() << endl;
        //}



        ND1_hasD = false;
        ND1_hasA = false;
        NE2_hasD = false;
        NE2_hasA = false;
        //found_something = false;
        //foundwhat.clear();
        //foundinres.clear();
        //foundwithdist.clear();



      }

    }

  }


  return seq;
}

void writeResSeq(std::ostream &os, std::vector<std::string> seq) {
  
  for (unsigned int i = 0; i < seq.size(); i++) {
    os << setw(6);
    
    if (i % 10 == 0 && i > 0) {
      os << endl;
    }
    os << seq[i];
  }
  os << endl;
  
}

void writePDB(gio::InPDB &myPDB, std::vector<std::string> seq){

  // should be included in OutPDB.cc ....
  
  // we need to know the terminal residues:
  // Finding where to put end groups
  // By default, before the first residue and after the last residue
  
  vector<unsigned int> endposition;
  
  for(unsigned int i = 0; i < myPDB.numAtoms()-1; ++i) {
    if (myPDB.getChain(i) != myPDB.getChain(i + 1)) 
      endposition.push_back(myPDB.getResNumber(i));
  }
  endposition.push_back(myPDB.getResNumber(myPDB.numAtoms()-1));
  int counter = 0;

  FILE * pFile;
  pFile = fopen ("corrected.pdb","w");
  
  for(unsigned int i=0; i<myPDB.numAtoms(); ++i){
    
    string type = myPDB.getType(i);
    unsigned int serial;
    string atom;
    string resName;
    string chainID;
    unsigned int seqNo;
    char iCode;
    double x, y, z;
    double occupancy;
    double tempFactor;
    string element;
    string charge;
    
    if (myPDB.hasSerials()) {
      serial = myPDB.getSerial(i);
    }
    if(myPDB.hasAtomNames()) {
      atom = myPDB.getAtomName(i);
    }
    if(myPDB.hasResNames()) {
      resName = myPDB.getResName(i);
    }
    if(myPDB.hasChainIDs()) {
      chainID = myPDB.getChain(i);
    }
    if(myPDB.hasResSeqs()) {
      seqNo = myPDB.getResNumber(i);
    }
    if(myPDB.hasiCodes()) {
      iCode = myPDB.getICode(i);
    }
    if(myPDB.hasX() && myPDB.hasY() && myPDB.hasZ()) {
      gmath::Vec pos = myPDB.getAtomPos(i);
      x = pos[0];
      y = pos[1];
      z = pos[2];
    }
    if(myPDB.hasOccupancies()) {
      occupancy = myPDB.getOcc(i);
    }
    if(myPDB.hasTempFactors()) {
      tempFactor = myPDB.getB(i);
    }
    if(myPDB.hasElements()) {
      element = myPDB.getElement(i);
    }
    if(myPDB.hasCharges()) {
      charge = myPDB.getCharge(i);
    }
    
    unsigned int atomnum = myPDB.getSerial(i);
    string atomname = myPDB.getAtomName(i);
    
    //string resname = myPDB.getResName(i);
    string resname = seq[myPDB.getResNumber(i)-1];
    if(atomname == "CD1" && resname == "ILE") {
      atomname = "CD";
    }

    if (seqNo == endposition[counter] && atomname == "O") {
      atomname = "O1";
      counter++;
    }

    fprintf(pFile, "%-6s", type.c_str());
    if (myPDB.hasSerials()) {
      fprintf(pFile, "%5d", atomnum);
    } else {
      fprintf(pFile, "%5s", "");
    }
    fprintf(pFile, "  ");
    if (myPDB.hasAtomNames()) {
      fprintf(pFile, "%-4s", atomname.c_str());
    } else {
      fprintf(pFile, "%-4s", "");
    }
    if (myPDB.hasResNames()) {
      fprintf(pFile, "%-4s", resname.c_str());
    } else {
      fprintf(pFile, "%-4s", "");
    }
    if (myPDB.hasChainIDs()) {
      fprintf(pFile, "%1s", chainID.c_str());
    } else {
      fprintf(pFile, "%1s", "");
    }
    if (myPDB.hasResSeqs()) {
      fprintf(pFile, "%4d", seqNo);
    } else {
      fprintf(pFile, "%4s", "");
    }
    if (myPDB.hasiCodes()) {
      fprintf(pFile, "%1s", &iCode);
    } else {
      fprintf(pFile, "%1s", "");
    }
    fprintf(pFile, "   ");
    if (myPDB.hasX()) {
      fprintf(pFile, "%8.3f", x);
    } else {
      fprintf(pFile, "%8s", "");
    }
    if (myPDB.hasY()) {
      fprintf(pFile, "%8.3f", y);
    } else {
      fprintf(pFile, "%8s", "");
    }
    if (myPDB.hasZ()) {
      fprintf(pFile, "%8.3f", z);
    } else {
      fprintf(pFile, "%8s", "");
    }
    if (myPDB.hasOccupancies()) {
      fprintf(pFile, "%6.2f", occupancy);
    } else {
      fprintf(pFile, "%6s", "");
    }
    if (myPDB.hasTempFactors()) {
      fprintf(pFile, "%6.2f", tempFactor);
    } else {
      fprintf(pFile, "%6s", "");
    }
    fprintf(pFile, "%10s", "");
    if (myPDB.hasElements()) {
      fprintf(pFile, "%-2s", element.c_str());
    } else {
      fprintf(pFile, "%-2s", "");
    }
    if (myPDB.hasCharges()) {
      fprintf(pFile, "%-2s\n", charge.c_str());
    } else {
      fprintf(pFile, "%-2s\n", "");
    }

    
  }
  fclose(pFile);
}

void checkAdaptPDB(gio::InPDB &myPDB, utils::gromosAminoAcidLibrary &gaal) {
  
  // to avoid double WARNINGS for the same type of amino acid, lets remember
  // what previously has been found
  static set<string> notFoundResNames;
  
  map<std::string, gromosAminoAcid> aa = gaal.getAminoAcids();
  for(unsigned int a = 0; a < myPDB.numAtoms(); ++a) {
    string resName = myPDB.getResName(a);
    if (aa.find(resName) == aa.end()) {
      if (resName.size() > 3) {
        int start = resName.size() - 3;
        string newResName = resName.substr(start, 3);
        if (notFoundResNames.find(resName) == notFoundResNames.end()) {
          cerr << "WARNING: cannot find residue " << resName << " (pdb name) in the"
                  " amino acid library and use " << newResName << " instead.\n";
          notFoundResNames.insert(resName);
        }
        resName = newResName;
        if(aa.find(resName) == aa.end()) {
          stringstream msg;
          msg << "cannot find residue " << resName << " (pdb name) in the amino acid library";
          throw gromos::Exception("pdb2seq", msg.str());
        }
        myPDB.setResName(a, resName);
      }
    }
  }
}

