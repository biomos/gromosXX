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
 * @file gca.cc
 * Generate coordinates from (dihedral) angles
 */

/**
 * @page programs Program Documentation
 *
 * @anchor gca
 * @section gca Generate coordinates from (dihedral) angles
 * @author @ref co
 * @date 22-9-06
 *
 * Sometimes, one may want to modify a specifed molecular configuration such as
 * to obtain specified values of bond lengths, bond angles or dihedral angles. 
 * Program gca allows the user to do this. In addition, series of 
 * configurations can be generated in which the molecular properties of choice
 * are modified stepwise. If more than one property to change has been 
 * specified, configurations for all combinations of values will be generated, 
 * allowing for a systematic search of the property space. In order to fulfill 
 * the requested property values, program gca will
 * - for a bond length between atoms i and j, shift all atoms connected
 *   to j and onwards (default: \@mobile last) or i and backwards (\@mobile first)
 * - for a bond angle defined by atoms i,j,k, rotate all atoms connected to 
 *   k and onwards (default: \@mobile last) or i and backwards (\@mobile first) 
 *   around the axis through atom j and perpendicular to the i,j,k-plane;
 * - for a dihedral angle defined by atoms i,j,k,l, rotate all atoms 
 *   connected to k and l and onwards (default: \@mobile last) or i and j and 
 *   backwards (\@mobile first) around the axis through atoms j and k.
 * 
 * This procedure may lead to distortions elsewhere in the molecule if the atom
 * count is not roughly linear along the molecular structure, or if the 
 * specified properties are part of a cyclic structure. The program does not
 * check for steric clashes resulting from the modifications. The properties to
 * be modified are specified through a @ref PropertySpecifier, followed by 
 * either one additional argument (single value to be specified) or three
 * additional arguments (to generate a range of values).
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref PropertySpecifier "properties" to change&gt; </td></tr>
 * <tr><td> [\@mobile</td><td>&lt;@ref args::OutformatParser "output coordinates format"&gt;] </td></tr>
 * <tr><td> [\@outformat</td><td>&lt;which part of the molecule should be mobile: first or last (default)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;input coordinate file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  gca
    @topo       ex.top
    @pbc        r
    @prop       t%1:3,4,5,6%100%0%360
    @outformat  cnf
    @mobile     first
    @traj       exref.coo
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/OutCoordinates.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/Neighbours.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p, string mobile);
int in_property(Property &p, int i);
int in_atoms(AtomSpecifier &as, int i);
void check_existing_property(System &sys, Property &p);
int findBond(System &sys, utils::Property &pp);
int findAngle(System &sys, utils::Property &pp);
int findDihedral(System &sys, utils::Property &pp);
void move_atoms(System &sys, AtomSpecifier &as, Vec v);
void rotate_atoms(System &sys, AtomSpecifier &as, gmath::Matrix rot, Vec v);

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "prop" << "outformat" << "traj" << "mobile";

  string usage = "# "+ string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t@prop       <properties to change>\n";
  usage += "\t[@mobile       <which part of the molecule should be mobile: first or last (default)>]\n";
  usage += "\t[@outformat <output coordinates format>]\n";
  usage += "\t@traj       <input coordinate trajectory>\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());
    
    // parse boundary conditions
    Boundary *pbc= BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // input 
    InG96 ic(args["traj"]);
    
    //read output format
    string ext;
    OutCoordinates *oc = OutformatParser::parse(args, ext);
    oc->open(cout);  
    oc->select("ALL");
    
    // prepare the title
    ostringstream stitle;
    stitle << "gca has modified coordinates in " <<args["traj"]
	   << " such that";
    
    // which part of the system should be mobile
    string mobile = "last";
    if (args.count("mobile") > 0) {
      mobile=args["mobile"];
      if (mobile != "first" && mobile != "last") 
        throw gromos::Exception("gca",
              "@mobile argument '" + mobile + "' unknown.\n");
    }
    
    // what properties. 
    PropertyContainer props(sys, pbc);
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
	to=args.upper_bound("prop");
      if(iter==to)
	throw gromos::Exception("gca", 
				"no property specified");
      for(; iter!=to; ++iter)
	props.addSpecifier(iter->second.c_str());
      
    }
    vector<double> min, max, step;
    for(unsigned int i=0; i< props.size(); i++){
      double stepsize=0, minimum=0, maximum=0;
      if(props[i]->args().size()==1){
        // single value given
        stepsize=1;
	minimum=props[i]->args()[0].scalar();
	maximum=props[i]->args()[0].scalar();

	stitle << endl << props[i]->toTitle() << "\tis set to " << minimum;
      }
      else if(props[i]->args().size() ==3){
        // we'll do a series
        stepsize = props[i]->args()[0].scalar();
        minimum = props[i]->args()[1].scalar();
        maximum = props[i]->args()[2].scalar();
	
        stitle << endl << props[i]->toTitle()  << "\tgoes from " << std::setw(8) << minimum
	       << " to " << std::setw(8) << maximum << " with a step size of " << std::setw(8) << stepsize;
      }
      else
	throw gromos::Exception("gca",
				"properties: specify single value or step%min%max values");

      step.push_back(stepsize);
      min.push_back(minimum);
      max.push_back(maximum);
    }
    
    oc->writeTitle(stitle.str());

    // generate all combinations
    vector<vector<double> > combination;
    vector<double> single=min;
    unsigned int index=0;
    bool flag=true;
    
    while(flag && single[props.size()-1]<=max[props.size()-1]){
      combination.push_back(single);
      single[index]+=step[index];
      while(single[index]>max[index]){
        single[index]=min[index];
        index++;
        if(index>=props.size()){ flag=false; break; }
        single[index]+=step[index];
        if(single[index] <= max[index]) index=0;
      }
    }

/*
    for(unsigned int i=0; i< combination.size(); ++i){
      for(unsigned int j=0; j< combination[i].size(); ++j){
	cout << combination[i][j] << " ";
      }
      cout << endl;
    }
*/
  
    int stepnum = 0;
    double time = 0;
    
    // loop over all trajectories
    for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	
	(*pbc.*gathmethod)();
	// loop over the combinations
	for(unsigned int c=0; c<combination.size(); ++c){
	  
	  // loop over the properties
	  for(unsigned int i=0; i< props.size(); i++){
	    
	    // check whether this is actually a property we know from the topology
	    // For the program, this is not necessary, but if you define arbitrary
	    // properties, the decision on which atoms will move also gets 
	    // arbitrary
	    check_existing_property(sys, *props[i]);
	    
	    // store the target value
	    double target=combination[c][i];
	    
	    // calculate what the current value is      
	    double value=props[i]->calc().scalar();
	    
	    
	    // determine which atoms to move
	    AtomSpecifier as(sys);
	    atoms_to_change(sys, as, *props[i], mobile);
	    // int m=props[i]->atoms().mol(0);
	    
	    // now handle the three different cases
	    string type = props[i]->type();
	    if(type=="Distance"){
	      Vec v = ( props[i]->atoms().pos(1)
			- props[i]->atoms().pos(0)).normalize();
	      v *= (target-value);
	      move_atoms(sys, as, v);
	    }
	    else if(type=="Angle"){
	      Vec v1 = props[i]->atoms().pos(0)
		- props[i]->atoms().pos(1);
	      Vec v2 = props[i]->atoms().pos(2)
		- props[i]->atoms().pos(1);
	      Vec v3 = v1.cross(v2);
	      gmath::Matrix rot=fit::PositionUtils::rotateAround(v3, target-value);
	      rotate_atoms(sys, as, rot, props[i]->atoms().pos(1));
	    }
	    else if(type=="Torsion"){
	      Vec v = props[i]->atoms().pos(2)
		- props[i]->atoms().pos(1);
	      
	      gmath::Matrix rot=fit::PositionUtils::rotateAround(v, target-value);
	      rotate_atoms(sys, as, rot, props[i]->atoms().pos(2));
	    }
            else {
              ostringstream msg;
              msg << "cannot modify property of type " << type;
              throw gromos::Exception(argv[0], msg.str());
            }
	  }

	  oc->writeTimestep(stepnum, time);
	  *oc << sys;
	  
	  ++stepnum;
	  time += 1;
	  
	}
      }
      
    }
    
    oc->close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

int in_property(Property &p, int i)
{
  // checks whether the atom i, is part of the property definition
  for(unsigned int j=0; j< p.atoms().size(); j++){
    if(i==p.atoms().atom(j)) return 1;
  }
  return 0;
}

int in_atoms(AtomSpecifier &as, int i)
{
  // checks whether atom i appears in the atom specifier
  // As our properties are always real bonds, angles or dihedrals, we do
  // not care about the molecule number here.
  // The AtomSpecifier should maybe get a function like this on its own.
  for(unsigned int j=0; j< as.size(); j++){
    if(i==as.atom(j)) return 1;
  }
  return 0;
}

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p, string mobile)
{
  // Determines which atoms of the whole molecule are influenced by changing
  // the value of the property
  
  // for bonds and angles we start with the last atom
  // for torsions we start with the third atom

  // all atoms bounded to this atom but not part of the property 
  // will be modified

  // we make either the part before (mobile="first") or after (mobile="last")
  // the property mobile
  
  int end2, endt;
  set<int> a, new_a;
  if (mobile == "last") {
    end2=p.atoms().atom(p.atoms().size()-1);
    if(p.type()=="Torsion") endt=p.atoms().atom(2);
  } else if (mobile == "first") {
    end2=p.atoms().atom(0);
    if(p.type()=="Torsion") endt=p.atoms().atom(1);
  }
  a.insert(end2);
  if(p.type()=="Torsion") a.insert(endt);
  
  int m=p.atoms().mol(0);
  
  while(a.begin()!=a.end()){
    set<int>::const_iterator b=a.begin(), e=a.end();
    new_a.clear();
    for(; b!=e; ++b){
      as.addAtom(m, *b);
      Neighbours n(sys, m, *b);
      for(unsigned int i=0; i<n.size(); i++)
	if(!in_property(p, n[i]) && !in_atoms(as, n[i]))
	   new_a.insert(n[i]);
    }
    a.clear();
    a=new_a;
  }
}

void move_atoms(System &sys, AtomSpecifier &as, Vec v)
{
  // Move the atoms in the atom-specifier by the vector v
  int m, a;
  
  for(unsigned int i=0; i<as.size(); i++){
    m=as.mol(i);
    a=as.atom(i);
    sys.mol(m).pos(a) += v;
  }
}
void rotate_atoms(System &sys, AtomSpecifier &as, gmath::Matrix rot, Vec v)
{
  // Rotate the atoms in the atom-specifyer according to the Matrix rot.
  // before rotation the atoms are moved so that v is the origin, after
  // rotation the atoms are moved back by v again.
  int m, a;
  Vec t;
  for(unsigned int i=0; i<as.size(); i++){
    m = as.mol(i);
    a = as.atom(i);
    t=sys.mol(m).pos(a) - v;
    t = rot*t;
    sys.mol(m).pos(a) = t + v;
  }
}

void check_existing_property(System &sys, Property &p)
{
  // Checks whether the property p is defined in the topology.
  string type=p.type();
  if(type=="Distance"){
    findBond(sys, p);
  }
  else if(type=="Angle"){
    findAngle(sys, p);
  }
  else if(type=="Torsion"){
    findDihedral(sys, p);
  }
  else{
    
    throw(gromos::Exception("gca", 
			    "Unknown property type"+p.type()));
  }
}

int findBond(System &sys, utils::Property &pp){
  // searches for the bond, defined by the property in the topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b,f=0;
  if(pp.atoms().mol(0)==pp.atoms().mol(1))
      m=pp.atoms().mol(0);
  else
    throw gromos::Exception("gca",
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms().atom(0)<pp.atoms().atom(1)){
    a=pp.atoms().atom(0);
    b=pp.atoms().atom(1);
  }
  else{
    a=pp.atoms().atom(1);
    b=pp.atoms().atom(0);
  }
  BondIterator bi(sys.mol(m).topology());
  while(bi && f==0)
    if(bi()[0]==a && bi()[1]==b) return 1;
    else ++bi;
  throw gromos::Exception("gca", 
			  " Bond not found in topology: "+pp.toTitle());
}

 int findAngle(System &sys, utils::Property &pp){
  // searches for the angle, defined by the property in the topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b,c,f=0;
  if(pp.atoms().mol(0)==pp.atoms().mol(1)&&pp.atoms().mol(0)==pp.atoms().mol(2))
      m=pp.atoms().mol(0);
  else
    throw gromos::Exception("gca",
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms().atom(0)<pp.atoms().atom(2)){
    a=pp.atoms().atom(0);
    b=pp.atoms().atom(1);
    c=pp.atoms().atom(2);
  }
  else{
    a=pp.atoms().atom(2);
    b=pp.atoms().atom(1);
    c=pp.atoms().atom(0);
  }
  AngleIterator ai(sys.mol(m).topology());
  while(ai&&f==0)
      if(ai()[0]==a && ai()[1]==b && ai()[2]==c) return 1;
    else ++ai;
  throw gromos::Exception("gca", 
        " Angle not found in topology: "+pp.toTitle());
}
    
int findDihedral(System &sys, utils::Property &pp){
  // searches for the (improper) dihedral, defined by the property in the 
  // topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b, c, d, f=0;
  if(pp.atoms().mol(0)==pp.atoms().mol(1) &&
     pp.atoms().mol(0)==pp.atoms().mol(2) &&
     pp.atoms().mol(0)==pp.atoms().mol(3))
      m=pp.atoms().mol(0);
  else
    throw gromos::Exception("gca", 
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms().atom(1)<pp.atoms().atom(2)){
    a=pp.atoms().atom(0);
    b=pp.atoms().atom(1);
    c=pp.atoms().atom(2);
    d=pp.atoms().atom(3);
  }
  else{
    a=pp.atoms().atom(3);
    b=pp.atoms().atom(2);
    c=pp.atoms().atom(1);
    d=pp.atoms().atom(0);
  }
  DihedralIterator di(sys.mol(m).topology());
  while(di&&f==0)
    if(di()[0]==a && di()[1]==b && di()[2]==c && di()[3]==d) return 1;
    else ++di;
  //Maybe we have an improper
  ImproperIterator ii(sys.mol(m).topology());
  while(ii&&f==0) 
    if(ii()[0]==a && ii()[1]==b && ii()[2]==c && ii()[3]==d) return 1;
    else ++ii;
  throw gromos::Exception("gca", 
			  " (improper) Dihedral not found in topology: "+pp.toTitle());
}

