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
 * @file ion.cc
 * replace water molecules by ions
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ion
 * @section ion replace water molecules by ions
 * @author @ref co
 * @date 11-6-07
 *
 * When simulating a charged solute in solution, one may wish to include
 * counter-ions in the molecular system in order to obtain a neutral system, or
 * a system with a specific ionic strength. The program ion can replace solvent
 * molecules by atomic ions by placing the
 * ion at the position of the first atom of a solvent molecule. Substitution of
 * solvent molecules by positive or negative ions can be performed by selecting
 * the solvent positions with the lowest or highest Coulomb potential, respectively,
 * or by random selection. In order to prevent two ions being placed too
 * close together, a sphere around each inserted ion can be specified from which
 * no solvent molecules will be substituted by additional ions. In addition, the user can
 * specify specific water molecules that should not be considered for
 * replacement.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> [\@positive</td><td>&lt;number&gt; &lt;ionname&gt; &lt;residue name (optional)&gt;] </td></tr>
 * <tr><td> [\@negative</td><td>&lt;number&gt; &lt;ionname&gt; &lt;residue name (optional)&gt;] </td></tr>
 * <tr><td> [\@potential</td><td>&lt;cutoff for potential calculation&gt;] </td></tr>
 * <tr><td> [\@random</td><td>&lt;random seed&gt;] </td></tr>
 * <tr><td> [\@exclude</td><td>&lt;solvent @ref AtomSpecifier "molecules" to be excluded&gt;] </td></tr>
 * <tr><td> [\@mindist</td><td>&lt;minimum distance between ions&gt;] </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ion
    @topo         ex.top
    @pbc          r
    @positive     10 NA+ NA
    @negative     10 CL- CL
    @potential    1.4
    @exclude      s:1-33
    @mindist      0.35
    @pos          exref.coo
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/OutG96S.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Energy.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

vector<int> selectPotential(utils::AtomSpecifier &sol, 
			    unsigned int num[2], double mindist,
			    utils::Energy &en, bound::Boundary *pbc);

vector<int> selectRandom(utils::AtomSpecifier &sol, 
			 int num, double mindist, 
			 int seed, bound::Boundary *pbc);

  
int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "positive" << "negative" << "potential"
         << "random" << "exclude" << "mindist" << "pos";

  string usage = argv[0];
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> <gather method>\n";
  usage += "\t[@positive  <number> <ionname> <residue name (optional)>]\n";
  usage += "\t[@negative  <number> <ionname> <residue name (optional)>]\n";
  usage += "\t[@potential <cutoff for potential calculation>]\n";
  usage += "\t[@random    <random seed>]\n";
  usage += "\t[@exclude   <solvent molecules to be excluded>]\n";
  usage += "\t[@mindist   <minimum distance between ions>]\n";
  usage += "\t@pos        <input coordinate file>\n";
 
try{
  Arguments args(argc, argv, knowns, usage);

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // read in coordinates
  InG96 ic(args["pos"]);
  ic.select("ALL");
  ic >> sys;
  
  // determine what to do.
  const int positive=0;
  const int negative=1;
  
  unsigned int num_ions[2]={0,0};
  std::string ion_names[2], res_names[2];
  {
    Arguments::const_iterator iter=args.lower_bound("positive"), 
      to = args.upper_bound("positive");
    if(iter!=to){
      num_ions[positive]=atoi(iter->second.c_str());
      ++iter;
    }
    if(iter!=to){
      ion_names[positive] = iter->second;
      ++iter;
    }
    else if(num_ions[positive])
      throw gromos::Exception("ion", "No name specified for positive ions");
    if(iter!=to)
      res_names[positive] = iter->second;
    else
      res_names[positive] = ion_names[positive];
    iter=args.lower_bound("negative");
    to=args.upper_bound("negative");
    if(iter!=to){
      num_ions[negative]=atoi(iter->second.c_str());
      ++iter;
    }
    if(iter!=to){
      ion_names[negative] = iter->second;
      ++iter;
    }
    else if(num_ions[negative])
      throw gromos::Exception("ion", "No name specified for negative ions");
    if(iter!=to)
      res_names[negative] = iter->second;
    else
      res_names[negative] = ion_names[negative];
  }

  // random or potential-based positions?
  bool random = false, potential = false;
  int seed = 0;
  
  if(args.count("random")>0){
    random = true;
    seed = atoi(args["random"].c_str());
  }
  double cutoff = 0.0;
  if(args.count("potential")>0){
    potential = true;
    cutoff = args.getValue<double>("potential");
  }
  if(!random && !potential)
    throw gromos::Exception("ion", "Don't know what to do, please specify \"random\" or \"potential\" flag");
  if(random && potential)
    throw gromos::Exception("ion", "Please specify either \"random\" or \"potential\" flag");
  
  // read the minimum ion-ion distance
  double mindist = args.getValue<double>("mindist", false, 0.0);
  
  // and the excluded solvent atoms
  utils::AtomSpecifier exclude(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("exclude"), 
      to=args.upper_bound("exclude");
    for(;iter!=to; ++iter)
      exclude.addSpecifier(iter->second);
  }
  
  // create a new Molecule topologies for the ions
  gcore::MoleculeTopology mt[2];
  for(int t=positive; t<=negative; t++){
    if(num_ions[t]){
      gcore::AtomTopology at;
      at.setName(ion_names[t]);
      mt[t].addAtom(at);
      mt[t].setResNum(0,0);
      mt[t].setResName(0,res_names[t]);
    }
  }

  // create a list of all first atoms of the solvent
  utils::AtomSpecifier sol(sys);
  sol.addSolventType(sys.sol(0).topology().atom(0).name());

  // remove the excluded solvents
  for(unsigned int i=0; i<exclude.size(); i++)
    sol.removeAtom(exclude.mol(i), exclude.atom(i));
  
  if(sol.size()<num_ions[positive]+num_ions[negative])
    throw gromos::Exception("ion", "Less solvent molecules allowed for replacement than ions requested");
  
  // generate a vector of atoms to be removed
  vector<int> remove;
  
  if(potential){
    utils::Energy en(sys, gff, *pbc);
    en.setAtoms(sol);
    en.setCutOff(cutoff);
    remove = selectPotential(sol, num_ions, mindist, en, pbc);
  }
  
  if(random)
    remove = selectRandom(sol, num_ions[positive] + num_ions[negative],
			  mindist, seed, pbc);

  // now also create a set with the atoms to remove
  set<int> remove_set;
  for(unsigned int i=0; i<remove.size(); i++) 
    remove_set.insert(sol.atom(remove[i]));
  
  // create a new system. This will be some copying
  System nsys(it.system());
  //first copy the solute
  for(int m=0; m<sys.numMolecules(); m++){
    nsys.mol(m).initPos();
    for(int a=0; a< sys.mol(m).numAtoms(); a++)
      nsys.mol(m).pos(a) = sys.mol(m).pos(a);
  }
  
  // now, add the ions
  for(int t=positive; t<=negative; t++){
    for(unsigned int i=0; i<num_ions[t]; i++){
      Molecule mol(mt[t]);
      // initialize the coordinates
      mol.initPos();
      mol.pos(0) = *sol.coord(remove[i+t*num_ions[0]]);
      nsys.addMolecule(mol);
    }
  }

  // and copy the coordinates of the remaining solvent
  int nsa=sys.sol(0).topology().numAtoms();
  for(int i=0; i<sys.sol(0).numPos(); i+=nsa)
    if(remove_set.count(i)==0)
      for(int j=0; j<nsa; j++)
	nsys.sol(0).addPos(sys.sol(0).pos(i+j));
  //and finally the box
  nsys.box() = sys.box();
  
  
  // define an output coordinate
  OutG96S oc(cout);
  oc.select("ALL");
  ostringstream os;
  os << "ion has replaced " << num_ions[positive] +  num_ions[negative] 
     << " solvent molecules in " << args["pos"] << " by" << endl;
  if(num_ions[positive]) {
    os << num_ions[positive] << " positive ions (" << ion_names[positive]
       << ")";
    if(num_ions[negative]) os << " and ";
    else os << endl;
  }
  if(num_ions[negative])
    os << num_ions[negative] << " negative ions (" << ion_names[negative] 
       << ")" << endl;
  os << "Solvent molecules were selected ";
  if(random) os << "randomly";
  if(potential) os << "according to the highest / lowest" << endl
		   << "electrostatic potential";
  if(exclude.size()){
    os << endl;
    os << "Solvent atoms ";
    for(Arguments::const_iterator iter=args.lower_bound("exclude");
	iter != args.upper_bound("exclude"); ++iter)
      os << iter->second << " ";
    os << "were not considered";
  }
  
  oc.writeTitle(os.str());
  oc << nsys;
  oc.close();
}


  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





vector<int> selectPotential(utils::AtomSpecifier &sol, 
			 unsigned int num[2], double mindist, 
			 utils::Energy &en, bound::Boundary *pbc)
{
  vector<int> selected;
  // calculate the non-bonded interactions
  en.calcNb();
  
  gmath::Vec v1, v2;
  set<int> not_allowed;
  double charge = sol.charge(0);
  double mindist2 = mindist*mindist;
  int s[2] = { 1,-1 };
  
  for(int t=0; t<2; t++){
    double sign=charge*s[t];
    for(unsigned int i=0; i<num[t]; i++){
      double min=1e6;
      int min_index=-1;
    
      for(unsigned int j=0; j<sol.size(); j++){
	if(en.el(j) / sign < min && not_allowed.count(j)==0){
	  min=en.el(j)/sign;
	  min_index=j;
	}
      }
      if(min_index==-1)
	throw gromos::Exception("ion", "not enough solvents with low potential found");
      
      selected.push_back(min_index);

      v1 = sol.pos(min_index);
      // now put all solvents within mindist into not_allowed.
      for(unsigned int j=0; j<sol.size(); j++){
	v2 = pbc->nearestImage(sol.pos(j), v1, sol.sys()->box());
	if((v2 - sol.pos(j)).abs2() < mindist2){
	  not_allowed.insert(j);
	}
      }
    }
  }
  return selected;
}

vector<int> selectRandom(utils::AtomSpecifier &sol, 
			 int num, double mindist, 
			 int seed, bound::Boundary *pbc)
{
  vector<int> selected;
  vector<int> allowed;
  double mindist2=mindist*mindist;
  gmath::Vec v1, v2;
  
  for(unsigned int i=0; i<sol.size(); i++) allowed.push_back(i);
  
  srand(seed);
  for(int i=0; i<num; i++){
    double r=rand();
    int index =  int((r*allowed.size())/RAND_MAX);
    int remove = allowed[index];
    
    selected.push_back(remove);
    gmath::Vec v1=*sol.coord(remove);
    // now remove all solvents within mindist from allowed.
    for(unsigned int j=0; j<sol.size(); j++){
      v2 = pbc->nearestImage(*sol.coord(j), v1, sol.sys()->box());
      if((v2 - *sol.coord(j)).abs2() < mindist2)
	allowed.erase(allowed.begin() + index);
    }
  }
  return selected;
}



