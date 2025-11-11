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
 * @file crossdih_top.cc
 * Patches a topology in order to introduce cross dihedral terms
 */

/**
 * @page contrib Contrib program documentation
 *
 * @anchor crossdih_top
 * @section crossdih_top Introduce cross dihedral terms
 * @author @ref ns
 * @date 9-10-2009
 *
 * Takes a cross dihedral correction file and applies them to the topology.
 *
 * The correction file format:
 * @verbatim
TITLE
correction file
END
CROSSCORRECTIONS
# number of correction sets
2
# number of residue names affected followed by residue names
29 ALA ARG ARGN ASN ASN1 ASP ASPH CYS CYSH CYS1 CYS2 GLN GLU GLUH HISA HISB HISH HIS1 ILE LEU LYS LYSH MET PHE SER THR TRP TYP VAL
# dihedral type of phi and psi
39 C N CA C 40 N CA C N
# number of corrections
11
# the corrections
phi  1.2   1      0
phi  0.5   1   -120
phi  0.8   2   -180
phi  1.0   3      0
psi  0.8   1    150
psi  0.5   3    -90
psi  0.5   4     90
pp   0.5   1     90
pp   1.1   2      0
pp   0.5   3   -135
pp   0.5   4    120
# second set for glycine
1 GLY
39 C N CA C 40 N CA C N
6
phi   0.3   2      0
phi   0.8   3      0
psi   2.8   1   -180
psi   4.0   2   -180
psi   0.5   3   -180
pp    0.5   2      0
END
@endverbatim
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@corr</td><td>&lt;cross dihedral correction file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  com_top
    @topo    ex.top
    @corr    cross.dat
 @endverbatim
 *
 * <hr>
 */


#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <math.h>
#include <set>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

enum corrtype { corrphi, corrpsi, corrcross };

struct dihedral_corr {
  corrtype type;
  DihedralType new_type;
};

struct backbone_dih {
  int type;
  string i, j, k, l;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "corr";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@corr <cross dihedral correction file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    Ginstream corrfile(args["corr"]);

    vector<set<string> > residue_names;
    vector<backbone_dih> phi, psi;
    vector<vector<dihedral_corr> > corrections;

    int max_dihedral_type = 0;
    GromosForceField gff(it.forceField());

    for(int i = 0; i < gff.numDihedralTypes(); ++i) {
      max_dihedral_type = max<int>(max_dihedral_type, gff.dihedralType(i).code());
    }

    vector<string> block;
    while(corrfile.getblock(block)) {
      if (block[0] == "CROSSCORRECTIONS") {
        string b;
        gio::concatenate(block.begin()+1, block.end()-1, b);
        istringstream is(b);
        unsigned int numsets;
        is >> numsets;
        if (is.fail()) {
          throw gromos::Exception(argv[0], "Cannot read number of correction sets");
        }
        residue_names.resize(numsets);
        corrections.resize(numsets);
        phi.resize(numsets);
        psi.resize(numsets);
        for(unsigned int set = 0; set < numsets; ++set) {
          unsigned int numres;
          is >> numres;
          if (is.fail()) {
            throw gromos::Exception(argv[0], "Cannot read number of affected residues");
          }
          for (unsigned int i = 0; i < numres; ++i) {
            string resname;
            is >> resname;
            if (is.fail()) {
              throw gromos::Exception(argv[0], "Cannot read reside name");
            }
            residue_names[set].insert(resname);
          }
          is >> phi[set].type >> phi[set].i >> phi[set].j >> phi[set].k >> phi[set].l;
          is >> psi[set].type >> psi[set].i >> psi[set].j >> psi[set].k >> psi[set].l;
          if (is.fail()) {
            throw gromos::Exception(argv[0], "Cannot read phi and psi dihedral type");
          }
          --phi[set].type; --psi[set].type;
          unsigned int numcorr;
          is >> numcorr;
          if (is.fail()) {
            throw gromos::Exception(argv[0], "Cannot read number of corrections");
          }
          corrections[set].resize(numcorr);
          for(unsigned int i = 0; i < numcorr; ++i) {
            string type;
            is >> type;
            if (is.fail() || (type != "phi" && type != "psi" && type != "pp")) {
              throw gromos::Exception(argv[0], "Cannot read type");
            }
            if (type == "phi")
              corrections[set][i].type = corrphi;
            else if (type == "psi")
              corrections[set][i].type = corrpsi;
            else
              corrections[set][i].type = corrcross;

            double fc, pdl;
            int np;
            is >> fc >> np >> pdl;
            if (is.fail()) {
              throw gromos::Exception(argv[0], "Cannot read k, m or phase");
            }

            corrections[set][i].new_type = DihedralType(++max_dihedral_type, fc, cos(pdl * M_PI / 180.0), pdl, np);
          }
        } // for sets
      }
    }

    if (phi.empty())
      throw gromos::Exception(argv[0], "No corrections to apply.");

    std::ostringstream title;
    title << "cross dihedral corrected topology" << endl;
    title << "Applying " << phi.size() << " sets of corrections." << endl;

    LinearTopology lt(sys);

    set<Dihedral> dih_remove, dih_add;
    set<CrossDihedral> crossdih_add;
    for (unsigned int s = 0; s < corrections.size(); ++s) {
      vector<dihedral_corr> & corr = corrections[s];

      // loop over residues
      for (unsigned int resnum = 0; resnum < lt.resNames().size(); ++resnum) {
        const Dihedral *dih_phi = NULL, *dih_psi = NULL;

        // skip unaltered residues
        const string resname = lt.resNames()[resnum];
        if (!residue_names[s].count(resname)) continue;

        // find phi and psi
        for (set<Dihedral>::const_iterator it = lt.dihedrals().begin(), to = lt.dihedrals().end(); it != to; ++it) {
          const Dihedral & dih = *it;
          // skip dihedrals of different residues
          if (lt.resMap()[dih[1]] != int(resnum)) continue;

          if (dih.type() == phi[s].type &&
                  lt.atoms()[dih[0]].name() == phi[s].i &&
                  lt.atoms()[dih[1]].name() == phi[s].j &&
                  lt.atoms()[dih[2]].name() == phi[s].k &&
                  lt.atoms()[dih[3]].name() == phi[s].l) {
            dih_phi = &dih;
          }
          if (dih.type() == psi[s].type &&
                  lt.atoms()[dih[0]].name() == psi[s].i &&
                  lt.atoms()[dih[1]].name() == psi[s].j &&
                  lt.atoms()[dih[2]].name() == psi[s].k &&
                  lt.atoms()[dih[3]].name() == psi[s].l) {
            dih_psi = &dih;
          }
          if (dih_phi != NULL && dih_psi != NULL) break;
        } // for dihedrals

        if (dih_phi == NULL) {
          ostringstream msg;
          msg << "Unable to find phi dihedral in residue " << resnum+1 << " (" << resname << ")";
          cerr << msg.str() << std::endl;
          continue;
        }
        if (dih_psi == NULL) {
          ostringstream msg;
          msg << "Unable to find psi dihedral in residue " << resnum+1 << " (" << resname << ")";
          cerr << msg.str() << std::endl;
          continue;
        }

        dih_remove.insert(*dih_phi);
        dih_remove.insert(*dih_psi);
        for(vector<dihedral_corr>::const_iterator cit = corr.begin(), cto = corr.end(); cit != cto; ++cit) {
          if (cit->type == corrphi) {
            Dihedral newdih((*dih_phi)[0],(*dih_phi)[1],(*dih_phi)[2],(*dih_phi)[3]);
            newdih.setType(cit->new_type.code());
            dih_add.insert(newdih);
          } else if (cit->type == corrpsi) {
            Dihedral newdih((*dih_psi)[0],(*dih_psi)[1],(*dih_psi)[2],(*dih_psi)[3]);
            newdih.setType(cit->new_type.code());
            dih_add.insert(newdih);
          } else {
            CrossDihedral newdih((*dih_phi)[0],(*dih_phi)[1],(*dih_phi)[2],(*dih_phi)[3],(*dih_psi)[0],(*dih_psi)[1],(*dih_psi)[2],(*dih_psi)[3]);
            newdih.setType(cit->new_type.code());
            crossdih_add.insert(newdih);
          }
        } // apply corrections
      } // loop over residues
    } // loop over corrections

    // remove dihedrals
    for(set<Dihedral>::const_iterator it = dih_remove.begin(), to = dih_remove.end(); it != to; ++it) {
      lt.dihedrals().erase(lt.dihedrals().find(*it));
    }
    // add new dihedrals
    for(set<Dihedral>::const_iterator it = dih_add.begin(), to = dih_add.end(); it != to; ++it) {
      lt.dihedrals().insert(*it);
    }
    // add cross dihedrals
    for(set<CrossDihedral>::const_iterator it = crossdih_add.begin(), to = crossdih_add.end(); it != to; ++it) {
      lt.crossdihedrals().insert(*it);
    }
    // add new dihedral types
    for (unsigned int set = 0; set < corrections.size(); ++set) {
      vector<dihedral_corr> & corr = corrections[set];
      for(vector<dihedral_corr>::const_iterator cit = corr.begin(), cto = corr.end(); cit != cto; ++cit) {
        gff.addDihedralType(cit->new_type);
      }
    }

    System corrected_sys;
    lt.parse(corrected_sys);
    for(int i = 0; i < sys.numSolvents(); ++i)
      corrected_sys.addSolvent(sys.sol(i));

    // add the pressure and temperature groups
    for(int i = 0; i < it.system().numTemperatureGroups(); ++i) {
      corrected_sys.addTemperatureGroup(it.system().temperatureGroup(i));
    }
    for(int i = 0; i < it.system().numPressureGroups(); ++i) {
      corrected_sys.addPressureGroup(it.system().pressureGroup(i));
    }

    OutTopology ot(cout);

    ot.setTitle(title.str());
    ot.write(corrected_sys, gff);
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





