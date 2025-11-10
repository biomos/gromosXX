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
 * @file link_top.cc
 * Create a link in a topology file according to a linking building block
 */

/**
 * @page programs Program Documentation
 *
 * @anchor link_top
 * @section link_top Create a link in a topology file according to a linking building block
 * @author @ref co
 * @date 19-10-15
 *
 * For branched systems, it may be very cumbersome to create crosslinks in a topology
 * This program allows the user to apply a pre-defined link to the topology. The link 
 * is defined in a special building block file, which contains (after a TITLE block), 
 * MTBUILDBLLINK blocks. This block has the following layout:  
 *
 * @verbatim
MTBUILDBLLINK
# RNME
XYZ
# number of atoms
 7
#ATOM RES ANM  IACM    MASS        CGMICGM MAE MSAE
    1   1 CA     14  13.018    0.00000   1   2    2   6 
    2   1 CB     15  14.027    0.16000   1   2    5   6
    3   1 OG      0       0    0.00000   1   0
    4   1 HG      0       0    0.00000   1   0
    5   2 CB     15  14.027    0.16000   0   1    6 
    6   2 OG      4  15.994   -0.32000   1   0 
    7   2 HG      0       0    0.00000   1   0
# bonds
# NB
   2
#  IB   JB  MCB
    1    2    1
    2    6   12
# bond angles
# NBA
    2
#  IB   JB   KB  MCB
    1    2    6   12
    2    6    5   12
# improper dihedrals
# NIDA
    0
#  IB   JB   KB   LB  MCB
# dihedrals
# NDA
    3
#  IB   JB   KB   LB  MCB
    0    1    2    6    1
    1    2    6    5    1
    2    6    5    0    1
END
@endverbatim
 *
 * The atoms section of the building block contains all atoms that are involved in the 
 * link. The second column specifies that these atoms are to be found in the first or
 * second residue of the link. The atoms are identified in the original topology by the
 * residue sequence number indicated in the input (\@linking) and the name of the atom
 * according to the MTBUILDBLLINK.
 *
 * In a first step, link_top, removes all atoms for which the IAC is 0. All references to
 * these atoms (exclusions, bonds, angles, etc.) are removed from the topology.
 * Next, the remaining atoms in the link defition get updated: the values in the topology
 * of IAC, mass, charge, and charge group get replaced by whatever is indicated in the 
 * MTBUILDBLLINK block. The exclusions of the original topology (without the removed 
 * atoms) remain, and the exclusions that are specified in the MTBUILDBLLINK block are 
 * added.
 *
 * Covalent interactions that need to be changed are also specified in the MTBUILDBLLINK
 * block. The program only allows the user to specify bonds, angles and improper dihedral
 * angles that are referring to atoms that are all part of the link specification. Any 
 * bonds, angles and improper dihedral angles that were present in the topology for these
 * atoms will be removed and the newly defined interactions are added to the topology. 
 * For dihedral angles, the program allows the user to refer to the first and/or last 
 * atom to be represented by a number 0. For these atoms, the program will search in the 
 * topology for an atom that is bound to the second or third atom, respectively and
 * assign the dihedral angle to this atom. Any dihedral angles that were already defined
 * for this group is replaced. Multiple dihedral angles for the same set of four atoms may
 * be added.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology&gt; </td></tr>
 * <tr><td> \@links</td><td>&lt;file with the MTBUILDBLLINK blocks&gt; </td></tr>
 * <tr><td> \@linking</td><td>&lt;residue 1&gt; &lt;residue 2&gt; &lt;link name&gt; </td></tr>
 * </table>
 *

 * Example:
 * @verbatim
  link_top
    @topo   ex.top
    @links  link.mtb
    @linking  2  4  XYZ
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/BbLink.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InLinkingBlock.h"
#include "../src/gio/OutTopology.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo" << "links" << "linking";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@links  <linking block file>\n";
  usage += "\t@linking <residue 1> <residue 2> <link name>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    // read in the linking block file
    InLinkingBlock ilb(args["links"]);

    std::vector<gcore::BbLink> links(ilb.links());

    // read in what we need to do
    std::vector<std::vector<int> > res;
    std::vector<std::string> name;
    if(args.count("linking") % 3 != 0)
      throw gromos::Exception("link_top",
                              "Error parsing input. @linking flag expects multiple of three arguments");

    for(Arguments::const_iterator iter=args.lower_bound("linking"),
          to=args.upper_bound("linking"); iter!=to;++iter) {
      int idum=atoi(iter->second.c_str());
      std::vector<int> tmp;
      tmp.push_back(idum-1);
      iter++;
      idum = atoi(iter->second.c_str());
      tmp.push_back(idum-1);
      res.push_back(tmp);
      iter++;
      name.push_back(iter->second);
    }

    // create a linear topology
    gcore::LinearTopology lt(sys);
    
    // loop over the tasks
    for(unsigned int i=0; i< name.size(); i++){
      //std::cerr << "performing the link " << res[i][0]+1 << " to " << res[i][1]+1 << " according to " << name[i] << std::endl;
      // find the appropriate link
      bool found=false;
      gcore::BbLink link;
      for(unsigned int j=0; j< links.size(); j++){
       	if(links[j].resName() == name[i]){
          found=true;
          link=links[j];
        }
      }
      if(!found){
        throw gromos::Exception("link_top",
                                "Link " + name[i] + " not found in the Linking file");
      }

      // find and flag the atoms that need to be removed
      for(int j=0; j < link.numAtoms(); j++){
        if(link.atom(j).iac() < 0){
          // we will need resiue res[i][j]  
          //std::cerr << "searching for residue " << res[i][link.linkRes(j)] +1 << " and atom " << link.atom(j).name() << std::endl;
          for(unsigned int k=0; k < lt.atoms().size(); k++){
            if(lt.resMap()[k] == res[i][link.linkRes(j)] && lt.atoms()[k].name() == link.atom(j).name()){
              lt.atoms()[k].setIac(-1);
            }
          }
        }
      }  
      // lt.removeAtoms();
   
      // now we make a map of the numbers in the LinkingBlock to the current
      // atom numbers in the reduced topology
      std::map<int, int> atommap;
      for(unsigned int j=0; j < (unsigned int)link.numAtoms(); j++){
        if(link.atom(j).iac() >= 0){
          bool found=false;
	  for(unsigned int k=0; k < lt.atoms().size(); k++){
            if(lt.resMap()[k] == res[i][link.linkRes(j)] && lt.atoms()[k].name() == link.atom(j).name()){
              //std::cerr << "found atom " << link.linkRes(j) << " " << link.atom(j).name() << std::endl;
              atommap[j]=k;
              found=true;
            }
          }
          if(!found){
            stringstream os;
            os << "Atom " << j+1 << " " << link.linkRes(j)+1 << " " << link.atom(j).name() << " in linking block not found in topology\n";
            throw gromos::Exception("link_top", os.str());
          }
        }
      }

      //for(std::map<int,int>::iterator iter=atommap.begin(); iter!=atommap.end(); ++iter){
      //  int j=iter->first;
        //std::cerr << "atommap " << j + 1 << " " << link.linkRes(j)+1  << " " << link.atom(j).name() << " is atom " << atommap[j]+1 << " " << lt.atoms()[atommap[j]].name() << std::endl;
      //}

      // Start applying changes
      //atoms
      for(std::map<int,int>::iterator iter=atommap.begin(); iter!=atommap.end(); ++iter){
        int j=iter->first;
        lt.atoms()[atommap[j]].setIac(link.atom(j).iac());
        lt.atoms()[atommap[j]].setMass(link.atom(j).mass());
        lt.atoms()[atommap[j]].setCharge(link.atom(j).charge());
        lt.atoms()[atommap[j]].setChargeGroup(link.atom(j).chargeGroup());
        // all existing exclusions remain, but the ones in the linking group are added
        gcore::Exclusion e=lt.atoms()[atommap[j]].exclusion();
	for(int k=0; k< link.atom(j).exclusion().size(); k++){
          // here's a risk: the exclusion may refer to something that is not in the atommap
          if(atommap.count(link.atom(j).exclusion().atom(k))==0){
            stringstream os;
            os << "Atom " << j+1 << " " << link.linkRes(j)+1 << " " << link.atom(j).name() 
               << " in linking block " << name[i] 
               << " has an exclusion that goes beyond the atoms in the linking block\n";
            throw gromos::Exception("link_top", os.str());
          } 
          e.insert(atommap[link.atom(j).exclusion().atom(k)]);
        }
        lt.atoms()[atommap[j]].setExclusion(e);
      }
      //bonds
      for(BondIterator bi(link); bi; ++bi){
        if(atommap.count(bi()[0])==0 || atommap.count(bi()[1])==0){
          stringstream os;	
          os << "Bond " << bi()[0]+1 << " " << bi()[1]+1  
             << " in linking block " << name[i] 
             << " refers to atoms that are not in the linking block, or no longer part of the topology\n";
          throw gromos::Exception("link_top", os.str());
        } 
        //std::cerr << "bond " << atommap[bi()[0]]+1 << " " << atommap[bi()[1]]+1 << " " << bi().type()+1 << std::endl;
        // see if the bond already exists in the topology, we remove it
        Bond b(0,1);
        bool found=false;
        for(std::set<Bond>::iterator iter=lt.bonds().begin(); iter!=lt.bonds().end(); ++iter){
          
          if((*iter)[0] == atommap[bi()[0]] && (*iter)[1] == atommap[bi()[1]]){
	    b=*iter;
            found=true;
            //std::cerr << "found bond" << std::endl;
          }
        }
        if(found)
          lt.bonds().erase(b);

        // and add the new one
        Bond bn(atommap[bi()[0]], atommap[bi()[1]]);
        bn.setType(bi().type());
        lt.bonds().insert(bn);
        
      }

      //angles
      for(AngleIterator ai(link); ai; ++ai){
        if(atommap.count(ai()[0])==0 || 
           atommap.count(ai()[1])==0 || 
           atommap.count(ai()[2])==0){
          stringstream os;	
          os << "Angle " << ai()[0]+1 << " " << ai()[1]+1 << " " << ai()[2]+1  
             << " in linking block " << name[i] 
             << " refers to atoms that are not in the linking block, or no longer part of the topology\n";
          throw gromos::Exception("link_top", os.str());
        } 
        // see if the angle already exists in the topology, we remove it
        Angle a(0,1,2);
        bool found=false;
        for(std::set<Angle>::iterator iter=lt.angles().begin(); iter!=lt.angles().end(); ++iter){
          if((*iter)[0] == atommap[ai()[0]] &&
             (*iter)[1] == atommap[ai()[1]] &&
             (*iter)[2] == atommap[ai()[2]]){
            a=*iter;
            found=true;
          }
        } 
        if(found)
          lt.angles().erase(a);

        // and add the new one
        Angle an(atommap[ai()[0]], atommap[ai()[1]], atommap[ai()[2]]);
        an.setType(ai().type());
	lt.angles().insert(an);
      
      }
      //impropers
      for(ImproperIterator ii(link); ii; ++ii){
        if(atommap.count(ii()[0])==0 || 
           atommap.count(ii()[1])==0 || 
           atommap.count(ii()[2])==0 ||
           atommap.count(ii()[3])==0){
          stringstream os;	
          os << "Improper " << ii()[0]+1 << " " << ii()[1]+1 << " " << ii()[2]+1 << " " << ii()[3]+1  
             << " in linking block " << name[i] 
             << " refers to atoms that are not in the linking block, or no longer part of the topology\n";
          throw gromos::Exception("link_top", os.str());
        } 
        // see if the angle already exists in the topology, we remove it
        Improper im(0,1,2,3);
        bool found=false;
        for(std::set<Improper>::iterator iter=lt.impropers().begin(); iter!=lt.impropers().end(); ++iter){
          if((*iter)[0] == atommap[ii()[0]] &&
             (*iter)[1] == atommap[ii()[1]] &&
             (*iter)[2] == atommap[ii()[2]] &&
             (*iter)[3] == atommap[ii()[3]]){
            im=*iter;
            found=true;
          }
        } 
        if(found)
          lt.impropers().erase(im);

        // and add the new one
        Improper in(atommap[ii()[0]], atommap[ii()[1]], atommap[ii()[2]], atommap[ii()[3]]);
        in.setType(ii().type());
	lt.impropers().insert(in);
      } 

      //dihedrals
      //the dihedrals give two complications
      // 1. there can be more than one for any set of four atoms
      // 2. they will refer to atoms outside the linking block

      std::set<Dihedral> added;
      for(DihedralIterator di(link); di; ++di){
        //we create a local atommap
        std::vector<int> localmap(4, -1);
        // we only allow the first or the last atom to be unclear
        if(atommap.count(di()[1])==0 || 
           atommap.count(di()[2])==0){
          stringstream os;	
          os << "Dihedral " << di()[0]+1 << " " << di()[1]+1 << " " << di()[2]+1 << " " << di()[3]+1  
             << " in linking block " << name[i] 
             << "\nThe two central atoms should be part of the linking block and still be part of the final topology\n";
          throw gromos::Exception("link_top", os.str());
        } 
        localmap[1]=atommap[di()[1]];
        localmap[2]=atommap[di()[2]];
        // the following should identify atoms at position 0 that are bound to position 1, and not equal to position 2
        if(atommap.count(di()[0])==0){
          for(std::set<Bond>::iterator it=lt.bonds().begin(); it!=lt.bonds().end(); ++it){
            if((*it)[0] == localmap[1] && (*it)[1]!=localmap[2]) localmap[0]=(*it)[1];
            if((*it)[1] == localmap[1] && (*it)[0]!=localmap[2]) localmap[0]=(*it)[0];
          }
        }
	else
          localmap[0] = atommap[di()[0]];
        // the following should identify atoms at position 3 that are bound to position 2, and not equal to position 1
        if(atommap.count(di()[3])==0){
          for(std::set<Bond>::iterator it=lt.bonds().begin(); it!=lt.bonds().end(); ++it){
            if((*it)[0] == localmap[2] && (*it)[1]!=localmap[1]) localmap[3]=(*it)[1];
            if((*it)[1] == localmap[2] && (*it)[0]!=localmap[1]) localmap[3]=(*it)[0];
          }
        }
	else
          localmap[3] = atommap[di()[3]];
        //std::cerr << "localmap " << localmap[0] +1 << " " << localmap[1] +1 << " " << localmap[2] +1 << " " << localmap[3]+1 << std::endl;
        for(unsigned int k=0; k< localmap.size(); k++){
          if(localmap[k] == -1){
            stringstream os;	
            os << "Dihedral " << di()[0]+1 << " " << di()[1]+1 << " " << di()[2]+1 << " " << di()[3]+1  
               << " in linking block " << name[i] 
               << "\nCould not find atom " << di()[k] << " in topology\n"; 
            throw gromos::Exception("link_top", os.str());

          }
        }

        // see if the dihedral already exists in the topology, we remove it
        std::vector<Dihedral> dihs;
        bool found=false;
        for(std::set<Dihedral>::iterator iter=lt.dihedrals().begin(); iter!=lt.dihedrals().end(); ++iter){
          if((*iter)[0] == localmap[0] &&
             (*iter)[1] == localmap[1] &&
             (*iter)[2] == localmap[2] &&
             (*iter)[3] == localmap[3]){
            dihs.push_back(*iter);
            found=true;
          }
        } 
        if(found){
          for(unsigned int k=0; k < dihs.size(); k++){
            // don't remove a dihedral we have added ourselves earlier
            if(added.count(dihs[k])==0){
              lt.dihedrals().erase(dihs[k]);	
            }
          }
        }
        // and add the new one
        Dihedral dn(localmap[0], localmap[1], localmap[2], localmap[3]);
        dn.setType(di().type());
	lt.dihedrals().insert(dn);
        added.insert(dn);
      } 
    }
    lt.removeAtoms();
    lt.get14s();
    
    System syo = lt.parse(); 
    syo.addSolvent(sys.sol(0));
    

    // set the temperature and pressure groups
    {
      int a = 0;
      for (int m = 0; m < syo.numMolecules(); ++m) {
        a += syo.mol(m).numAtoms();
        syo.addPressureGroup(a);
        syo.addTemperatureGroup(a);
      }
    }
    // and write out the new topology

    OutTopology ot(cout);
    ostringstream os;
    os << "Modified topology based on " << args["topo"] << endl;
    os << "Applied following crosslinks according to definitions in " 
       << args["links"] << endl;
    for(unsigned int i=0; i< name.size(); i++){
      os << "- residues " << res[i][0]+1 << " to " << res[i][1]+1 << " according to " 
         << name[i] << std::endl;
    }
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
   
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}



