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
 * @file gca_ener.cc
 * combines the functionality of @ref gca and @ref ener 
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor gca_ener
 * @section gca_ener combines the functionality of @ref gca and @ref ener
 * @author @ref co
 * @date 31-10-08
 *
 * Program gca_ener combines the functionalities of the programs @ref gca and 
 * @ref ener to avoid a disk-consuming step of writing an extra trajectory 
 * file. It allows the user to modify some bonded interactions as described 
 * for @ref gca and subsequently calculate the (non)bonded energies as
 * described for @ref ener.
 *
 * See the original programs for more detailed documentation.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref PropertySpecifier "properties" to change&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider&gt; </td></tr>
 * <tr><td> \@prop_ener</td><td>&lt;properties to calculate energies for&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field correction&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa for reaction field correction&gt; </td></tr>
 * <tr><td> \@RFex</td><td>&lt;calculate RF correction for excluded atoms: on/off&gt; </td></tr>
 * <tr><td> \@soft</td><td>&lt;soft @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;lam&gt; &lt;a_lj&gt; &lt;a_c&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input coordinate file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  gca_ener
    @topo        ex.top
    @pbc         r
    @prop        d%1:1,2%0.15
    @atoms       1:a
    @prop_ener   d%1:1,2
    @time        0 1
    @cut         1.4
    @eps         62.0
    @kap         0.0
    @RFex        on
    @soft        1:5
    @softpar     0.5 1.51 0.5
    @traj        ex.tr
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
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
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
#include "../src/utils/groTime.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Energy.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p);
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
  knowns << "topo" << "pbc" << "prop" << "traj" << "atoms" << "prop_ener" 
         << "time" << "cut" << "eps" << "kap" << "soft" << "softpar" << "RFex";

  string usage = "# "+ string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t@prop       <properties to change>\n";
  usage += "\t@atoms      <atoms to consider>\n";
  usage += "\t@prop_ener  <properties to calculate energies for>\n";
  usage += "\t@time       <time> <dt>\n";
  usage += "\t@cut        <cut-off distance>\n";
  usage += "\t@eps        <epsilon for reaction field correction>\n";
  usage += "\t@kap        <kappa for reaction field correction>\n";
  usage += "\t@RFex       <calculate RF for excluded atoms: on/off>\n";
  usage += "\t@soft       <soft atoms>\n";
  usage += "\t@softpar    <lam> <a_lj> <a_c>\n";
  usage += "\t@traj       <input coordinate file>\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc= BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // input 
    InG96 ic(args["traj"]);
    
    // prepare the title
    ostringstream stitle;
    stitle << "gca has modified coordinates in " <<args["traj"]
	   << " such that";
    
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
  
    // prepare the energy parts
    //   get simulation time
    Time time(args);
    // declare the energy class
    Energy en(sys, gff, *pbc);
 
    // RF for excluded atoms?
    {
      std::string s=args.getValue<string>("RFex",1,"on");
      if(s=="off")
	en.setRFexclusions(false);
      else
	en.setRFexclusions(true);
    }



    //  set atoms
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      Arguments::const_iterator to=args.upper_bound("atoms");
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	atoms.addSpecifier(spec);
      }
    }
    en.setAtoms(atoms);
    
    // set properties
    PropertyContainer prop_ener(sys, pbc);
    {
      Arguments::const_iterator iter=args.lower_bound("prop_ener");
      Arguments::const_iterator to=args.upper_bound("prop_ener");
      for(;iter!=to;iter++){
	string p=iter->second.c_str();
	prop_ener.addSpecifier(p);
      }
    }
    en.setProperties(prop_ener);

    // set non-bonded parameters
    //   get cut-off distance
    en.setCutOff(args.getValue<double>("cut", false, 1.4));

    en.setRF(args.getValue<double>("eps", false, 1.0),
            args.getValue<double>("kap", false, 0.0));
    // get soft atom list
    AtomSpecifier soft(sys);
    {
      bool lsoft = false;
      Arguments::const_iterator iter=args.lower_bound("soft");
      Arguments::const_iterator to=args.upper_bound("soft");
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	soft.addSpecifier(spec);
	lsoft = true;
      }
      //  get softnef parameters
      std::vector<double> softpar = args.getValues<double>("softpar", 3, lsoft,
              Arguments::Default<double>() << 0.0 << 0.0 << 0.0);
      if (lsoft)
        en.setSoft(soft, softpar[0], softpar[1], softpar[2]);
    }
    // print titles
    cout << "# Time"
	 << "           vanderwaals"
	 << "         electrostatic"
	 << "                 Total"
	 << endl;

    // declare some variables for averaging
    int num_frames=0;
    double vdw=0.0;
    double crf=0.0;
    double tot=0.0;

    int stepnum = 0;
    
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
	ic >> sys >> time;
	
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
	    atoms_to_change(sys, as, *props[i]);
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

	  // now we can calculate the energies
	  en.calc();
	  
	  // print any ouput you like
	  cout.precision(10);
	  cout.setf(ios::right, ios::adjustfield);
	  cout << setw(6) << time
	       << setw(22) << en.vdw()
	       << setw(22) << en.el()
	       << setw(22) << en.tot()
	       << endl;
	  
	  //store some averages
	  vdw+=en.vdw();
	  crf+=en.el();
	  tot+=en.tot();
	 
	  num_frames++;
	  
	  ++stepnum;
	  
	}
      }
      
    }
    // print out average energies
    if(num_frames>1){
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << endl << "# ave."
	   << setw(22) << vdw/num_frames 
	   << setw(22) << crf/num_frames
	   << setw(22) << tot/num_frames
	   << endl;
    }
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
  for(int j=0; j< p.atoms().size(); j++){
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
  for(int j=0; j< as.size(); j++){
    if(i==as.atom(j)) return 1;
  }
  return 0;
}

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p)
{
  // Determines which atoms of the whole molecule are influenced by changing
  // the value of the property

  // there are always two ends to a property. Regardless of the order of the
  // atoms we change the last end.
  int end2=p.atoms().atom(p.atoms().size()-1);

  // all atoms bounded to this atom but not part of the property 
  // will be modified
  set<int> a, new_a;
  a.insert(end2);
  // for bonds and angles we start with the last atom
  // for torsions we start with the third atom
  if(p.type()=="Torsion") a.insert(p.atoms().atom(2));
  else a.insert(end2);
  
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
  
  for(int i=0; i<as.size(); i++){
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
  for(int i=0; i<as.size(); i++){
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

