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
 * @file pairlist.cc
 * prints out a list of all particles within a cutoff distance from an atom
 * or a point in space
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pairlist
 * @section pairlist prints out a list of all particles within a cutoff 
 * @author @ref co
 * @date 25.8.2006
 *
 * Program pairlist determines all particles within user specified cutoffs from
 * a given reference point. The reference point can either be a single
 * @ref utils::AtomSpecifier "atom", a @ref utils::VectorSpecifier "vector" 
 * (preceded by the keyword "vector") or a set of three cartesian coordinates.
 * The output can be written in the same style as program @ref atominfo to 
 * allow usage as an @ref utils::AtomSpecifier "atomspecifier" itself.
 *
 * The program can produce two pairlists at the time, one shortrange and one
 * longrange. It will also print out a list of particles that occur in the
 * longrange pairlist only. The pairlist determination can be done on an atomic
 * basis or based  on charge groups.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr> 
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt;</td></tr>
 * <tr><td> \@coord</td><td>&lt;coordinaes to base the list on&gt; </td></tr>
 * <tr><td> \@refpos</td><td>&lt;reference @ref AtomSpecifier "atom"&gt; or &lt;vector&gt;</td></tr>
 * <tr><td> [\@cutp</td><td>&lt;small cutoff&gt;] </td></tr>
 * <tr><td> [\@cutl</td><td>&lt;large cutoff&gt;] </td></tr>
 * <tr><td> [\@type</td><td>&lt;ATOMIC (default) or CHARGEGROUP&gt;] </td></tr>
 * <tr><td> [\@atominfo</td><td>(write in atominfo style)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  atominfo
    @topo    ex.top
    @pbc     r
    @coord   exref.co
    @refpos  1:1
#   @refpos  vector atom(va(cog,1:2,4,5))
#   @refpos  0.1 0.5 1.5
    @cutp    0.8
    @cutl    1.4
    @type    CHARGEGROUP
    @atominfo
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

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Value.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace utils;

using namespace std;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "coord" << "cutl" << "cutp" << "refpos" 
         << "type" << "atominfo";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <boundary type> <gathermethod>\n";
  usage += "\t@coord     <coordinates to base the list on>\n";
  usage += "\t@refpos    <atoms or vector>\n";
  usage += "\t[@cutp     <small cutoff>]\n";
  usage += "\t[@cutl     <large cutoff>]\n";
  usage += "\t[@type     <ATOMIC (default) or CHARGEGROUP>]\n";
  usage += "\t[@atominfo <write in atominfo style>]\n";
    
  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // define in and output coordinates
    InG96 ic(args["coord"]);
    ic.select("ALL");
    
    ic >> sys;
    // just in case we use a virtual atom as the reference position: do gather
    (*pbc.*gathmethod)();
    
    // read in the cutoffs
    double cutp = args.getValue<double>("cutp", false, 0.0);
    double cutl = args.getValue<double>("cutl", false, 0.0);
    bool do_short=false, do_long=false, do_diff=false;
    if(cutp!=0.0) do_short=true;
    if(cutl!=0.0) do_long =true;
    do_diff = do_short && do_long;
    if(do_diff && cutp > cutl)
      throw gromos::Exception("pairlist", "cutp should be shorter than cutl");
    
    if(!do_short && !do_long)
      throw gromos::Exception("pairlist", 
			      "Please specify at least one of cutp or cutl");
    
    // read in the refpos
    gmath::Vec ref;
    AtomSpecifier ref_as(sys);
    VectorSpecifier vs(sys,pbc);
    
    bool ref_is_atom=false;
    
    int numMol=sys.numMolecules();

    if(args.count("refpos")==2){
      Arguments::const_iterator iter=args.lower_bound("refpos");
      if(iter->second != "vector")
	throw gromos::Exception("pairlist", 
				"Use the keyword 'vector' to use a vector specifier");
      iter++;
      VectorSpecifier vs(sys,pbc,iter->second);
      ref=vs();
    }
    	
    if(args.count("refpos")==3){
      Arguments::const_iterator iter=args.lower_bound("refpos");
      for(int i=0;iter != args.upper_bound("refpos"); ++iter, ++i)
	ref[i]=atof(iter->second.c_str());
    }
    if(args.count("refpos")==2 || args.count("refpos")==3){
      MoleculeTopology mt;
      AtomTopology at;
      at.setName("ref");
      at.setChargeGroup(1);
      mt.addAtom(at);
      sys.addMolecule(mt);
      sys.mol(numMol).initPos();
      sys.mol(numMol).pos(0) = ref;
      ref_as.addAtom(numMol,0);
    }
    if(args.count("refpos")==1){
      ref_as.addSpecifier(args["refpos"]);
      if(ref_as.size()!=1)
	throw gromos::Exception("pairlist",
		"only one atom should be specified with refpos");
      ref_is_atom=true;
    }

    // read in the type
    std::string t="ATOMIC";
    if(args.count("type")>0){
      if(args["type"]=="CHARGEGROUP") t=args["type"];
      else if(args["type"]!="ATOMIC") throw gromos::Exception("pairlist",
		"only ATOMIC and CHARGEGROUP are accepted as type");
    }
    
    SimplePairlist pp(sys, *pbc, cutp);
    SimplePairlist pl(sys, *pbc, cutl);
    AtomSpecifier  diff(sys);
    
    if(do_short){
      pp.setAtom(ref_as.mol(0), ref_as.atom(0));
      pp.setType(t);
      pp.calc();
      if(ref_is_atom)
	pp.addAtom(ref_as.mol(0), ref_as.atom(0));
    }
    if(do_long){
      pl.setAtom(ref_as.mol(0), ref_as.atom(0));
      pl.setType(t);
      pl.calc();
      if(ref_is_atom)
	pl.addAtom(ref_as.mol(0), ref_as.atom(0));
    }
    if(do_diff){
      for(unsigned int i=0; i<pl.size(); i++)
	if(pp.findAtom(pl.mol(i), pl.atom(i))==-1)
	  diff.addAtom(pl.mol(i), pl.atom(i));
    }

    // now produce the output
    cout << "# Making atom lists based on " << args["coord"] << endl << endl;
    cout << "# Reference point is " << (*ref_as.coord(0))[0] << " " 
	 <<  (*ref_as.coord(0))[1] << " " <<(*ref_as.coord(0))[2];
    if(ref_is_atom){
      cout << " (coordinates of atom ";
      if(ref_as.mol(0)<0) cout << "s"; else cout << ref_as.mol(0) + 1;
      cout << ":" << ref_as.atom(0) + 1 << ")";
    }
    cout << endl << endl;
    cout << "# Using a";
    if(t=="ATOMIC") cout << "n atomic";
    else cout << " charge group based";
    cout << " cutoff criterion" << endl << endl;
    
    if(do_short){
      cout << "# Within the short range ( r < " << cutp << " ) there are "
	   << pp.size() << " elements:" << endl << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tpairlist: short range\n";
	cout << "\nEND\n";
	cout << "ATOMS\n";
	
	cout << "#"
	     << setw(12) << "Atom"
	     << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	     << setw(10) << "Residue"
	     << setw(10) << "Atom"
	     << setw(12) << "Integer"
	     << setw(10) << "Charge" << endl;
	cout << "#"
	     << setw(12) << "Specifier"
	     << setw(10) << "number"
	     << setw(10) << "number"
	     << setw(10) << "name"
	     << setw(10) << "name"
	     << setw(12) << "Atom Code"
	     << endl;
	
	for(unsigned int i=0; i < pp.size(); ++i){
	  cout << setw(13) << pp.toString(i)
	       << setw(10) << pp.gromosAtom(i)+1
	       << setw(10) << pp.resnum(i)+1
	       << setw(10) << pp.resname(i)
	       << setw(10) << pp.name(i)
	       << setw(12) << pp.iac(i)+1
	       << setw(10) << pp.charge(i)
	       << endl;
	}
	cout << "END\n";
      }
      else{
	vector<string> s = pp.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
      }
      cout << endl << endl;
    }
    if(do_long){
      cout << "# Within the long range ( r < " << cutl << " ) there are "
	   << pl.size() << " elements:" << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tpairlist: long (and short) range\n";
	cout << "\nEND\n";
	cout << "ATOMS\n";
	
	cout << "#"
	     << setw(12) << "Atom"
	     << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	     << setw(10) << "Residue"
	     << setw(10) << "Atom"
	     << setw(12) << "Integer"
	     << setw(10) << "Charge" << endl;
	cout << "#"
	     << setw(12) << "Specifier"
	     << setw(10) << "number"
	     << setw(10) << "number"
	     << setw(10) << "name"
	     << setw(10) << "name"
	     << setw(12) << "Atom Code"
	     << endl;

	for(unsigned int i=0; i < pl.size(); ++i){
	  cout << setw(13) << pl.toString(i)
	       << setw(10) << pl.gromosAtom(i)+1
	       << setw(10) << pl.resnum(i)+1
	       << setw(10) << pl.resname(i)
	       << setw(10) << pl.name(i)
	       << setw(12) << pl.iac(i)+1
	       << setw(10) << pl.charge(i)
	       << endl;
	}
	cout << "END\n";
      }
      else{
	vector<string> s = pl.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
      }
      cout << endl << endl;
    }
    if(do_diff){
      cout << "# Within the shell ( " << cutp << " < r < " << cutl 
	   << " ) there are " << diff.size() << " elements:" << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tpairlist: long (only) range\n";
	cout << "\nEND\n";
	cout << "ATOMS\n";
	
	cout << "#"
	     << setw(12) << "Atom"
	     << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	     << setw(10) << "Residue"
	     << setw(10) << "Atom"
	     << setw(12) << "Integer"
	     << setw(10) << "Charge" << endl;
	cout << "#"
	     << setw(12) << "Specifier"
	     << setw(10) << "number"
	     << setw(10) << "number"
	     << setw(10) << "name"
	     << setw(10) << "name"
	     << setw(12) << "Atom Code"
	     << endl;


	for(unsigned int i=0; i < diff.size(); ++i){
	  cout << setw(13) << diff.toString(i)
	       << setw(10) << diff.gromosAtom(i)+1
	       << setw(10) << diff.resnum(i)+1
	       << setw(10) << diff.resname(i)
	       << setw(10) << diff.name(i)
	       << setw(12) << diff.iac(i)+1
	       << setw(10) << diff.charge(i)
	       << endl;
	}
	cout << "END\n";
	
      }
      else{
	vector<string> s = diff.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
      }
      cout << endl << endl;
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

