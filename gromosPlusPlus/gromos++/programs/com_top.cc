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
 * @file com_top.cc
 * Combine molecular topology files into one
 */

/**
 * @page programs Program Documentation
 *
 * @anchor com_top
 * @section com_top Combine molecular topology files into one
 * @author @ref co
 * @date 7-6-07
 *
 * To generate molecular topology files for the use in simulations of e.g. 
 * (macro)molecular complexes, or mixtures containing several solutes and/or 
 * (co)solvents, it is usually convenient to merge existing molecular topology
 * files. Program com top combines multiple topologies into one new topology.
 * 
 * The user has to specify which molecular topologies to be merged (use prefix
 * 'n:' before the file name to repeat one topology n times), and from  which 
 * file the force field parameters and the solvent have to be taken. The
 * resulting molecular topology file is written out to the standard output. 
 *
 * The program can also be used for topology file format conversion. The 
 * argument \@inG96 converts GROMOS96 topologies to the current format. On
 * the other hand \@outG96 converts topologys in the current format to the
 * GROMOS96 format.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology files&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;index number of molecular topology file to take parameters from&gt; </td></tr>
 * <tr><td> \@solv</td><td>&lt;index number of molecular topology file to take solvent from&gt; </td></tr>
 * </table>
 *
 *
 * Example 1:
 * @verbatim
  com_top
    @topo    ex.top 7:cl.top
    @param   1
    @solv    1
 @endverbatim
 *
 * Example 2 (format conversion):
 * @verbatim
  com_top
    @topo    ex.top
    @param   1
    @solv    1
    @inG96
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/VirtualAtom.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "param" << "solv";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology files>\n";
  usage += "\t@param <index number of molecular topology file to take parameters from>\n";
  usage += "\t@solv  <index number of molecular topology file to take solvent from>\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    
    // read and check the command line arguments
    //
    // first the topology
    int numtop = 1;
    if (args.count("topo") <= 0) {
      throw gromos::Exception("com_top", "needs at least one topology\n" + usage);
    } else {
      numtop = args.count("topo");
    }
    // then the @param and @solv flags
    int parnum = 1, solnum = 1;
    {
      string firsttopo = args.lower_bound("topo")->second;
      size_t pos = firsttopo.find(":");
      if (pos != string::npos) {
        firsttopo = firsttopo.substr(pos + 1);
      }
      if (args.count("param") > 0) {
        parnum = args.getValue<int>("param");
      } else {
        cerr << "NOTE: force field parameters taken from " << firsttopo << endl;
      }
      if (parnum < 1 || parnum > numtop) {
        if (numtop > 1) {
          stringstream msg;
          msg << "bad value for @param, read " << parnum << ", allowed is an integer from 1.." << numtop;
          throw gromos::Exception(argv[0], msg.str());
        } else {
          throw gromos::Exception(argv[0], "bad value for @param which must be 1");
        }
      }

      if (args.count("solv") > 0) {
        solnum = args.getValue<int>("solv");
      } else {
        cerr << "NOTE: solvent taken from " << firsttopo << endl;
      }
      if (solnum < 1 || solnum > numtop) {
        if (numtop > 1) {
          stringstream msg;
          msg << "bad value for @solv, read " << solnum << ", allowed is an integer from 1.." << numtop;
          throw gromos::Exception(argv[0], msg.str());
        } else {
          throw gromos::Exception(argv[0], "bad value for @solv which must be 1");
        }
      }
    }

    System sys;

    //ugly solution to the not yet implemented '=' for force fields
    string paramname, toponame, s;
    std::string::size_type s_it;
    int repeat = 1, topnum=0, oldtopnum=-1;
    
    ostringstream title;
    title << "COM_TOP: Combined topology using:\n";
    
   
    int totNumAt=0;    
    int totVirtAt=0;
    vector<int> virt_per_mol; 
    vector<int> tot_per_mol;
    for(Arguments::const_iterator iter=args.lower_bound("topo"),
         to=args.upper_bound("topo"); iter!=to; ++iter){
      oldtopnum=topnum+1;
      
      s=iter->second;
      s_it=s.find(':');
      
      if(s_it == string::npos){
	toponame=s;
	repeat=1;
      }
      else{
	toponame=s.substr(s_it+1,s.size());
	repeat=atoi(s.substr(0,s_it).c_str());
      
	if(repeat==0){
	  throw gromos::Exception("com_top", "using " + s.substr(0,s_it) + " in " + iter->second +
				  " is not allowed\n         use program con_top to exchange only force field parameters\n");
 }
      }
      
      topnum+=repeat;
      
      // read topology
      InTopology it(toponame);
      for(int i=0; i<repeat; i++){
        // Directly add pressure and temperature groups
        for (int j = 0; j < it.system().numTemperatureGroups(); j++) {
          sys.addTemperatureGroup(it.system().temperatureGroup(j) + totNumAt);
        }

        for (int j = 0; j < it.system().numPressureGroups(); j++) {
          sys.addPressureGroup(it.system().pressureGroup(j) + totNumAt);
        }

        // Add molecules 
        for (int j = 0; j < it.system().numMolecules(); j++) {
          sys.addMolecule(it.system().mol(j));
          cerr << "totVirtAt " << totVirtAt << endl;
          virt_per_mol.push_back(totVirtAt);
          tot_per_mol.push_back(totNumAt);
        }
        for (unsigned int j = 0; j < it.system().vas().numVirtualAtoms(); j++) {
          std::vector<int> conf;
          for(int k=0; k< it.system().vas().atom(j).conf().size(); k++){
            conf.push_back(it.system().vas().atom(j).conf().gromosAtom(k) + totNumAt);
          }
          sys.addVirtualAtom(conf, it.system().vas().atom(j).type(),
                                   0.1, 0.153, it.system().vas().iac(j), 
                                   it.system().vas().charge(j),
                                   it.system().vas().exclusion(j),
                                   it.system().vas().exclusion14(j)); 
          totVirtAt++;
          cerr << "added totVirtAt " << totVirtAt << endl;
        }
        for (int j = 0; j < it.system().numMolecules(); j++) {
          totNumAt += it.system().mol(j).numAtoms();
        } // molecules

      } // repeat

      if(solnum <= topnum && solnum >= oldtopnum)
	sys.addSolvent(it.system().sol(0));
      if(parnum <= topnum && parnum >= oldtopnum)
        paramname=toponame;
      if(topnum!=oldtopnum)
	title << setw(4) << oldtopnum << " .. " << setw(4) << topnum;
      else
	title << setw(12) << topnum;
      
      title << " : " << toponame << endl;
    }
    // if there are virtual atoms, there may be exclusions that need
    // adjustment
    for(int i=0; i < sys.numMolecules(); i++){
      for(int j=0; j < sys.mol(i).numAtoms(); j++){
        for(int k=0; k < sys.mol(i).topology().atom(j).exclusion().size(); k++){
          int current=sys.mol(i).topology().atom(j).exclusion().atom(k);
          if(current >= sys.mol(i).numAtoms()) {
            sys.mol(i).topology().atom(j).exclusion().erase(current);
            sys.mol(i).topology().atom(j).exclusion().insert(current + virt_per_mol[i] + totNumAt -tot_per_mol[i] - sys.mol(i).numAtoms());
          }
        } 
        for(int k=0; k < sys.mol(i).topology().atom(j).exclusion14().size(); k++){
          int current=sys.mol(i).topology().atom(j).exclusion14().atom(k);
          if(current >= sys.mol(i).numAtoms()) {
            sys.mol(i).topology().atom(j).exclusion14().erase(current);
            sys.mol(i).topology().atom(j).exclusion14().insert(current + virt_per_mol[i] + totNumAt -tot_per_mol[i] - sys.mol(i).numAtoms());
          }
        } 
      }
    }

    InTopology it(paramname);
    title << "Parameters from " << parnum 
          << ", solvent from " << solnum;
    
    OutTopology ot(cout);
   
    ot.setTitle(title.str());
    ot.write(sys,it.forceField());
  } catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




