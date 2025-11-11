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
 * @file prep_xray_le.cc
 * Creates input-file for xray local elevation
 */

/**
 * @page programs Program Documentation
 *
 * @anchor prep_xray_le
 * @section prep_xray_le Creates input-file for X-ray local elevation
 * @author @ref ns
 * @date 12-02-2010
 *
 * Program prep_xray_le creates a X-ray local elevation file. It takes the side chains
 * of the residues contained in @ref AtomSpecifier atoms. The sidechains are defined
 * in a special file (\@library). It should contain the following block:
 * @verbatim
LESIDECHAIN
  # name  dim  atom names
  ARG     4    N CA CB CG CA CB CG CD CB CG CD NE CG CD NE CZ
  ASN     2    N CA CB CG CA CB CG OD1
 END
@endverbatim
 * As the atom names define (dim) dihedral angles they have to be a multiple of four.
 * The local elevation parameters (force constant, the number of bins of the grid,
 * the functional form switch, the width of the potential and its cutoff) are specified
 * using \@leparam. The X-ray parameters (@f$R^0@f$ threshold and cutoff for
 * @f$R_\mathrm{real}@f$ calculation) are specified using \@xrayparam.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library that specifies dihedrals&gt; </td></tr>
 * <tr><td>[\@leparam</td><td>&lt;local elevation parameters&gt;]</td></tr>
 * <tr><td>[\@xrayparam</td><td>&lt;xray parameter&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  prep_xray_le
    @topo          ex.top
    @atoms         1:res(SER:1)
    @library       def.lib
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <ios>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace std;
using namespace utils;

int findAtomResidue(System & sys, int mol, int res, string name) {
  int atom_start = -1, atom_end = -1;
  for(int i = 0; i < sys.mol(mol).numAtoms(); ++i) {
    if (res == sys.mol(mol).topology().resNum(i)) {
      if (atom_start == -1) atom_start = i;
      atom_end = i;
    }
  }
  for (int i = atom_start; i < atom_end; ++i) {
    if (sys.mol(mol).topology().atom(i).name() == name) {
      return i - atom_start;
    }
  }
  return -1;
}

int main(int argc, char *argv[]) {
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@library        <specification file for dihedrals>\n";
  usage += "\t@atoms          <atoms to consider>\n";
  usage += "\t[@leparam       <force bins func width cutoff>]\n";
  usage += "\t[@xrayparam     <threshold cutoff>]\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "library" << "atoms" << "leparam" << "xrayparam";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    AtomSpecifier atoms(sys);
    for (Arguments::const_iterator it = args.lower_bound("atoms"),
            to = args.upper_bound("atoms"); it != to; ++it) {
      atoms.addSpecifier(it->second);
    }
    if (atoms.empty())
      throw gromos::Exception(argv[0], "No atoms specified.");
    atoms.sort();
    { // loop over the atoms and include the first atom of the residue
      set<pair<int,int> > res;
      for(unsigned int i = 0; i < atoms.size(); ++i)
        res.insert(pair<int,int>(atoms.mol(i), atoms.resnum(i)));
      atoms.clear();
      for(set<pair<int,int> >::const_iterator it = res.begin(), to = res.end();
              it != to; ++it) {
        ostringstream spec;
        spec << (it->first + 1) << ":res(" << (it->second + 1) << ":1)";
        atoms.addSpecifier(spec.str());
      }
      atoms.sort();
    }

    Ginstream libfile(args["library"]);
    vector<string> buffer;
    /*
    vector<int> backbone_spec;
    {
      libfile.getblock(buffer);
      if (buffer[0] != "LEBACKBONE")
        throw gromos::Exception(argv[0],
              "library file does not contain a LEBACKBONE block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception(argv[0], "library file file " + libfile.name() +
              " is corrupted. No END in LEBACKBONE block. Got\n" + buffer[buffer.size() - 1]);
      if (buffer.size() != 3)
        throw gromos::Exception(argv[0], "LEBACKBONE block has more than one line");
      
      istringstream line(buffer[1]);
      int dim;
      if (!(line >> dim))
        throw gromos::Exception(argv[0], "LEBACKBONE block: cannot read dimension.");
      backbone_spec.resize(dim);
      for(unsigned int i = 0; i < backbone_spec.size(); ++i)
        line >> backbone_spec[i];
      if (line.fail())
        throw gromos::Exception(argv[0], "LEBACKBONE block: cannot read all dihedral angle numbers.");
    }*/

    map<string, vector<string> > sidechain_spec;
    {
      libfile.getblock(buffer);
      if (buffer[0] != "LESIDECHAIN")
        throw gromos::Exception(argv[0],
              "library file does not contain a LESIDECHAIN block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception(argv[0], "library file file " + libfile.name() +
              " is corrupted. No END in LESIDECHAIN block. Got\n" + buffer[buffer.size() - 1]);

      for (vector<string>::const_iterator it = buffer.begin() + 1, to = buffer.end() - 1; it != to; ++it) {
        istringstream line(*it);
        string res;
        if (!(line >> res))
          throw gromos::Exception(argv[0], "LESIDECHAIN block: cannot read residue name.");
        int dim;
        if (!(line >> dim))
          throw gromos::Exception(argv[0], "LESIDECHAIN block: cannot read dimension.");
        sidechain_spec[res].resize(dim*4);
        for (unsigned int i = 0; i < sidechain_spec[res].size(); ++i)
          line >> sidechain_spec[res][i];
        if (line.fail())
          throw gromos::Exception(argv[0], "LESIDECHAIN block: cannot read all dihedral angle atom names for residue " + res + ".");
      }
    }

    double force_constant = 100.0;
    int bins = 36;
    int func = 0; // trunc poly
    double width = 1.0;
    double cutoff = 1.0;
    if (args.count("leparam") > 0) {
      stringstream ss;
      for(Arguments::const_iterator it = args.lower_bound("leparam"),
              to = args.upper_bound("leparam"); it != to; ++it)
        ss << it->second << " ";
      ss >> force_constant >> bins >> func >> width >> cutoff;
      if (ss.fail())
        throw gromos::Exception(argv[0], "invalid @leparam argument.");
    }

    vector<double> xrayparam = args.getValues<double>("xrayparam", 2, false,
          Arguments::Default<double>() << 0.2 << 0.25);
    double r_thres = xrayparam[0];
    double r_cut = xrayparam[1];

    // loop over residues within atomspec
    vector<string> res_desc;
    vector<vector<int> > res_indices;
    vector<set<int> > res_unique_indices;
    for (unsigned int r = 0; r < atoms.size(); ++r) {
      string resname = atoms.resname(r);
      // skip residue if not present in lib
      if (sidechain_spec.find(resname) == sidechain_spec.end())
        continue;
      
      ostringstream desc;
      desc << setw(5) << atoms.mol(r)+1 << setw(5) << atoms.resnum(r)+1 << setw(5) << resname;
      res_desc.push_back(desc.str());
      // found it. Loop over this residue
      vector<int> indices;
      for(unsigned int i = 0; i < sidechain_spec[resname].size(); ++i) {
        int index = findAtomResidue(sys, atoms.mol(r), atoms.resnum(r), sidechain_spec[resname][i]);
        if (index == -1) {
          ostringstream msg;
          msg << "Cannont find atom " << sidechain_spec[resname][i] << " in residue " << resname;
          throw gromos::Exception(argv[0], msg.str());
        }
        indices.push_back(atoms.gromosAtom(r) + index);
      }
      set<int> unique_indices;
      unique_indices.insert(indices.begin(), indices.end());
      
      res_indices.push_back(indices);
      res_unique_indices.push_back(unique_indices);
    }

    assert(res_desc.size() == res_indices.size());
    assert(res_desc.size() == res_unique_indices.size());

    // do umbrella weights first
    cout << "XRAYUMBRELLAWEIGHT\n" << "# UMB THRES CUTOFF  ATOMS\n";
    cout.precision(5);
    for(unsigned int i = 0; i < res_desc.size(); ++i) {
      cout << "#" << res_desc[i] << "\n";
      cout << setw(5) << i << setw(10) << r_thres << setw(10) << r_cut;
      for(set<int>::const_iterator it = res_unique_indices[i].begin(), to = res_unique_indices[i].end(); it != to; ++it)
        cout << setw(5) << (*it + 1);
      cout << "\n";
    }
    cout << "END\n";

    // and now the local elevation specification
    cout << "LOCALELEVSPEC\n";
    for(unsigned int r = 0; r < res_desc.size(); ++r) {
      cout << "#" << res_desc[r] << "\n";
      assert(res_indices[r].size() % 4 == 0);
      for(unsigned int i = 0; i < res_indices[r].size(); i += 4) {
        cout << setw(5) << r << setw(5) << res_indices[r][i] + 1
                << setw(5) << res_indices[r][i+1] + 1
                << setw(5) << res_indices[r][i+2] + 1
                << setw(5) << res_indices[r][i+3] + 1 << "\n";
      }
    }
    cout << "END\n";
    // and finally the umbrellas
    cout << "LEUSBIAS\n# NUMUMB\n"
            << res_desc.size() << "\n";
    cout.precision(8);
    for(unsigned int r = 0; r < res_desc.size(); ++r) {
      cout << "#" << res_desc[r] << "\n"
              << "# NLEPID NDIM CLES\n";
      unsigned int dim = res_indices[r].size() / 4;
      cout << setw(5) << r << setw(5) << dim << setw(16) << force_constant << "\n"
              << "# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)\n";
      for(unsigned int i = 0; i < dim; ++i) {
        cout << setw(5) << 1
                << setw(5) << func
                << setw(13) << width
                << setw(13) << cutoff
                << setw(5) << bins
                << setw(13) << 0.0
                << setw(13) << 0.0 << "\n";
      }
      cout << "# NCONLE\n0\n# NVISLE ICONF\n";

    }
    cout << "END\n";

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
