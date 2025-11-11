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
 * @file predict_noe.cc
 * Predict possible NOE pairs from distance averages
 */


/**
 * @page programs Program Documentation
 *
 * @anchor predict_noe
 * @section predict_noe Predict NOE pairs from distance averages
 * @author @ref ns 
 * @date 15-07-2011
 *
 * Program predict_noe is used to predict possible NOE pairs from distance averages.
 * The program calculates and averages all possible NOE pairs from a trajectory.
 * The nuclei are taken from the atom specifier provided if those are found in the
 * NOE library file as well.
 * The averaging is carried out as
 * @f[ \overline{r} = \left<r^{-p}\right>^{-\frac{1}{p}} @f],
 * where @f$ p @f$ can be either 1, 3 or 6.
 * Distances above a threshold level (\@filter) are discarded in the final output.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;Atoms to consider&gt; </td></tr>
 * <tr><td> \@lib</td><td>&lt;NOE specification library&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;coordinate trajectory files&gt; </td></tr>
 * <tr><td> [\@dish</td><td>&lt;carbon-hydrogen distance; default: 0.1 nm&gt;] </td></tr>
 * <tr><td> [\@disc</td><td>&lt;carbon-carbon distance; default: 0.153 nm&gt;] </td></tr>
 * <tr><td> [\@averaging</td><td>&lt;averaging power, 1, 3 or 6; default 6&gt;] </td></tr>
 * <tr><td> [\@filter</td><td>&lt;discard NOE's above a certain distance [nm]; default 10000 nm&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  predict_noe
    @topo          ex.top
    @atoms         1:a
    @pbc           r
    @lib           noelib.45a3
    @traj          ex.trc
    @dish          0.1
    @disc          0.153
    @averaging     6
    @filter        0.8
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/VirtualAtom.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Noe.h"
#include "../src/gcore/Molecule.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/VirtualAtoms.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;
using namespace gmath;

int main(int argc, char *argv[]) {

  // Usage string

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type>\n";
  usage += "\t@atoms        <atoms to consider>\n";
  usage += "\t@lib          <NOE specification library>\n";
  usage += "\t@traj         <trajectory files>\n";
  usage += "\t[@dish        <carbon-hydrogen distance; default: 0.1 nm>]\n";
  usage += "\t[@disc        <carbon-carbon distance; default: 0.153 nm>]\n";
  usage += "\t[@filter      <discard NOE's above a certain distance [nm]; default 10000 nm>]\n";
  usage += "\t[@averaging   <averaging 1, 3, 6>\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "pbc" << "filter" << "atoms" << "lib" << "dish" << "disc" << "averaging" << "traj";

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);

    //read in filter threshold
    const double filter = args.getValue<double>("filter", false, 10000.0);
    // Read in and create the NOE library
    Ginstream nff(args["lib"]);
    vector<string> buffer;
    nff.getblock(buffer);

    vector<Noelib> noelib;
    parse_noelib(buffer, noelib);
    nff.close();
    

    //try for disc and dish
    double dish = args.getValue<double>("dish", false, 0.1);
    double disc = args.getValue<double>("disc", false, 0.153);
    
    int averaging = args.getValue<int>("averaging", false, 6);
    if (averaging != 1 && averaging != 3 && averaging != 6) {
      throw gromos::Exception(argv[0], "averaging has to be 1, 3 or 6.");
    }
    
    AtomSpecifier specatoms(sys);

    //get rmsd atoms
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");

      for (; iter != to; iter++) {
        specatoms.addSpecifier(iter->second);
      }
    }
    if (specatoms.size() == 0)
      throw gromos::Exception(argv[0], "No atoms specified!");
    
    vector<VirtualAtom*> atoms;
    vector<pair<string, unsigned int> > atom_description;
    for(unsigned int i = 0; i < specatoms.size(); ++i) {
      vector<Noelib>::const_iterator it = noelib.begin(),
              to = noelib.end();
      for(; it != to; ++it) {
        if (specatoms.name(i) == it->gratomname && specatoms.resname(i) == it->resname) {
          ostringstream desc;
          desc << setw(5) << specatoms.mol(i) + 1 << " " << setw(5) << specatoms.atom(i) + 1
                  << " " << setw(5) << specatoms.name(i);
          // found
          vector<VirtualAtom*> v = getvirtual(specatoms.gromosAtom(i), it->NOETYPE,
                  it->NOESUBTYPE, sys, dish, disc);
          atoms.insert(atoms.end(), v.begin(), v.end());
          atom_description.push_back(pair<string, unsigned int>(desc.str(), v.size()));
          break;
        }
      }
    }
    const unsigned int atom_size = atoms.size();
    {
      vector<pair<string, unsigned int> >::const_iterator it = atom_description.begin(),
              to = atom_description.end();
      for (unsigned int i = 0; it != to; ++it) {
        for(unsigned int j = 0; j < it->second; ++j, ++i) {
          cout << "# " << setw(8) << i + 1 << ": " << it->first 
                  << " variant " << j + 1 << endl;
        }
      }
    }
    
    InG96 ic;
    unsigned int numFrames = 0;
    vector<vector<double> > matrix(atom_size, vector<double>(atom_size, 0.0));
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);


      // loop over all frames
      while (!ic.eof()) {

        numFrames++;
        ic >> sys;
        
#pragma omp parallel for
        for(unsigned int i = 0; i < atom_size; ++i) {
          const Vec & ri = atoms[i]->pos();
          for(unsigned int j = i + 1; j < atom_size; ++j) {
            const Vec & rj = atoms[j]->pos();
            const double dist2 = (ri - pbc->nearestImage(ri, rj, sys.box())).abs2();
            const double dist2i = 1.0 / dist2;
            switch(averaging) {
              case 1:
                matrix[i][j] += sqrt(dist2i);
                break;
              case 3:
                matrix[i][j] += sqrt(dist2i) * dist2i;
                break;
              case 6:
                matrix[i][j] += dist2i * dist2i * dist2i;
                break;
              default:
                throw gromos::Exception(argv[0], "averaging method not implemented.");
            }
          }
        }
        
      }
    }
    cout.precision(8);
    
    for(unsigned int i = 0; i < atom_size; ++i) {
      for(unsigned int j = i + 1; j < atom_size; ++j) {
        matrix[i][j] = pow(matrix[i][j] / double(numFrames), -1.0 / averaging);
        if (matrix[i][j] > filter)
          continue;
        
        cout << setw(8) << i + 1 << setw(8) << j + 1 << setw(15) << matrix[i][j]
                << endl;
      }
    }
    
    for(unsigned int i = 0; i < atom_size; ++i)
      delete atoms[i];
    
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
