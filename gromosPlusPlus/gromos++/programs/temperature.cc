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
 * @file temperature.cc
 * Calculate the temperature for different sets of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor temperature
 * @section temperature Calculates the temperature for different sets of atoms
 * @author nb
 * @date 1.6.2012
 *
 * Program temperature will calculate the temperature for different sets of atoms,
 * as specified by the atomspecifier(s). Multiple sets of atoms can be specified by
 * white-space separated Atomspecifiers. For each of the sets one dof value is
 * expected.
 * 
 * You can find the number of degree of freedoms for a temperature group in the
 * md++ output file under "DEGREES OF FREEDOM" -> "DOF"
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>]
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "Atoms sets" &gt; </td></tr>
 * <tr><td> \@dofs</td><td>&lt;degrees of freedom &gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;velocity trajectory files&gt; </td></tr>
 * </table>
 *
 * @sa @ref AtomSpecifier
 *   
 * Example:
@verbatim
 temperature
    @topo       ex.top
    @time       0 0.1
    @atoms      1:1-30 1:35-56
    @dofs       26 34
    @traj       ex.trv

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

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Temperature.h"

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace std;
using namespace gmath;

int main(int argc, char **argv){
  
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atoms" << "time" << "dofs";

  string usage = "# " + string(argv[0]);
  usage += "# You can find the number of degree of freedoms for a temperature group in the\n";
  usage += "# md++ output file under \"DEGREES OF FREEDOM\" -> \"DOF\".\n";
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atoms      <atoms to consider for the temperature calculations>\n";
  usage += "\t@dofs       <the number of degree of freedoms for these atoms>\n";
  usage += "\t@traj       <trajectory files>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // Check the number of arguments
    bool error = false;
    if (!(args.count("atoms") > 0)){
      cerr << "No atoms specified with '@atoms'!\n";
      error = true;
    }
    if (!(args.count("dofs") > 0)){
      cerr << "The number of degree of freedoms is not specified with '@dofs'!\n";
      error = true;
    }
    if (!(args.count("traj") > 0)){
      cerr << "No trajectory files specified with '@traj'!\n";
      error = true;
    }
    if (!(args.count("topo") > 0)){
      cerr << "No No topology specified with '@topo'!\n";
      error = true;
    }
    if (args.count("atoms") != args.count("dofs")){
      cerr << "There are not as many degree of freedom numbers as atom specifiers!\n";
      error = true;
    }
    if (error){
      cerr << "\nErrors during arguments checking!\n";
      exit(1);
    }
  
    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // Temperature
    bool has_solvent = false;
    cout << "# Time ";
    int temp_group = 1;
    vector<Temperature> temps;
    ostringstream os;
    os << "#      ";
    Arguments::iterator as_i = args.lower_bound("atoms"),
            as_to = args.upper_bound("atoms");
    Arguments::iterator dofs_i = args.lower_bound("dofs");
    
    for (; as_i != as_to; as_i++, dofs_i++){
      stringstream ss_dof;
      ss_dof << dofs_i->second;
      double dofs;
      ss_dof >> dofs;
      temps.push_back(Temperature(AtomSpecifier(sys, as_i->second), dofs));
//      temp_group++;
      cout << "<Temperature Group " << temp_group++<< ">";
      os << "<"<< as_i->second << ">  "; 
      if (as_i->second.find("s:") < as_i->second.size()){
        has_solvent = true;
      }
    }
    cout << endl;
    cout << os.str()<< std::endl;
    
    if (has_solvent){
      cout << "# Solvent specified!" << endl;
      cout << "# It might be better to turn the solvents into solutes in the topology!" << endl;
    }
    
    
    
    // Process the files and generate the temperature
    Arguments::iterator file_i = args.lower_bound("traj"),
            file_to = args.upper_bound("traj");
    for (; file_i != file_to; file_i++){
      InG96 ic;
      ic.open(file_i->second);
      while(!ic.eof()){
        ic.select("ALL");
        ic >> sys >> time;
        if (!sys.hasVel){
          cerr << "No Velocity read!\n";
          exit(1);
        }

        cout << time;
        vector<Temperature>::iterator temp_i = temps.begin(),
            temp_to = temps.end();
        for (; temp_i != temp_to; temp_i++){
          cout << " \t" << temp_i->temperature(sys);
        }
        cout << endl;
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
