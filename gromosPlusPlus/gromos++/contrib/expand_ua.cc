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
 * @file expand_ua.cc
 * Replace united atoms with real atoms in topology
 */

/**
 * @page contrib Contrib Documentation
 *
 * @anchor expand_ua
 * @section expand_ua Replace united atoms in topology with real atoms
 * @author @ref P. Poliak
 * @date 20. 6. 2020
 *
 * This program searches for united atoms in building block and replaces them
 * with bare carbon and corresponding number of hydrogens. Added hydrogens
 * are put in the same chargegroup and exlusions are inherited from the carbon
 * atom. For atoms excluding the former united atom also hydrogen atoms are
 * put into exclusion list. Charge of the united atom stays on the carbon.
 * Then also C-H bonds are added. X-C-H angles are added for all neighbouring
 * atoms based on found C-X bonds. H-C-H angles are added as well.
 * TODO:
 *  1. Improper dihedral if chiral carbon is created
 *  2. Read UA IACs and their internal topology from library file
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@file</td><td>&lt;molecular building block file&gt; </td></tr>
 * <tr><td> \@build</td><td>&lt;name of the building block to modify&gt; </td></tr>
 * <tr><td> [\@names</td><td>&lt;atom names to be expanded into real atoms (by default all united atoms will be expanded)&gt;]</td></tr>
 * <tr><td> [\@number</td><td>&lt;new building block name (default is old name)&gt;]</td></tr>
 * </table>
 *
 *
 * Example 1:
 * @verbatim
  add_atom
    @file    53a6.mtb
    @build   PHE
    @names   CA CB
    @new     EXPHE
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <algorithm>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/OutBuildingBlock.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "file" << "build" << "names" << "new";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@file <mtb-file>\n";
  usage += "\t@build <buildingblock name>\n";
  usage += "\t[@names <names of atoms to expand (default ALL)>]\n";
  usage += "\t[@new <new building block name (default is original name)>]\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    InBuildingBlock ibb(args["file"]);
    BuildingBlock mtb(ibb.building());

    int index = mtb.findBb(args["build"]);
    BbSolute obb(mtb.bb(index - 1));
    BbSolute nbb;

    vector<string> names;
    for (Arguments::const_iterator iter = args.lower_bound("names"),
      to = args.upper_bound("names"); iter != to; ++iter) {
	    names.push_back(iter->second);
    }

    string new_name;
    if (args.count("new") != -1) new_name = args["new"];
    else new_name = args["build"];

    // IACs of united atoms, later to be read from libfile
    int ua_iacs[3] = {13, 14, 15};

    // Find united atoms
    int a_it = 0;
    while (a_it++ < obb.numAtoms() - 1) {
      AtomTopology at(obb.atom(a_it));
      if (names.size() == 0 || find(names.begin(), names.end(), at.name()) != names.end()) {
        // also treat IAC 18 (aliphatic CH2 in ring) as IAC 15
        if (at.iac() == 17) at.setIac(14);
        if (find(ua_iacs, ua_iacs + sizeof(ua_iacs) / sizeof(int), at.iac()) != ua_iacs + sizeof(ua_iacs) / sizeof(int)) {
          // Parts of the following code are unscrupulously stolen from add_atom (thanks Chris)
          BbSolute nbb;
          int start = a_it+1;
          int number = at.iac() - 12;
          obb.atom(a_it).setIac(11); // Turn UA carbon to bare carbon
          obb.atom(a_it).setMass(11);
          // Add to the same chargegroup
          int cg = obb.atom(a_it).chargeGroup();
          obb.atom(a_it).setChargeGroup(0); // 1 will be given to the last H atom
          string n = obb.atom(a_it).name();
          string unique_id;
          if (n.size() > 1)
            unique_id = n.substr(1);
          //prepare an empty atom
          AtomTopology atn;
          atn.setIac(19); // hydrogen bound to carbon
          atn.setMass(0);
          atn.setCharge(0.0);
          atn.setChargeGroup(0); // this we will adapt later
          atn.setName("H");

          //copy obb to nbb, but stop at start
          //first the previous exclusions
          for (int i = 0; i < obb.numPexcl(); i++) {
            Exclusion e;
            for (int j = 0; j < obb.pexcl(i).size(); j++)
              if (obb.pexcl(i).atom(j) < start - 1) {
                e.insert(obb.pexcl(i).atom(j));
              } else if (obb.pexcl(i).atom(j) == start - 1) { // if UA excluded, than also new H excluded
                e.insert(obb.pexcl(i).atom(j));
                for (int k = 0; k < number; k++) {
                  e.insert(obb.pexcl(i).atom(j) + k + 1);
                }
              }
              else {
                e.insert(obb.pexcl(i).atom(j) + number);
              }

            nbb.addPexcl(e);
          }

          //and the name
          nbb.setResName(obb.resName());

          //now, the atoms bis und mit start, pay attention to the exclusions.
          for (int i = 0; i < start; i++) {
            AtomTopology at(obb.atom(i));
            Exclusion e;
            for (int j = 0; j < at.exclusion().size(); j++)
              if (at.exclusion().atom(j) < start - 1)
                e.insert(at.exclusion().atom(j));
              else if (at.exclusion().atom(j) == start - 1) { // if UA excluded, than also new H excluded
                e.insert(at.exclusion().atom(j));
                for (int k = 0; k < number; k++) {
                  e.insert(at.exclusion().atom(j) + k + 1);
                }
              }
              else
                e.insert(at.exclusion().atom(j) + number);
            at.setExclusion(e);
            nbb.addAtom(at);
          }
          //add number atoms
          for (int i = 0; i < number; i++) {
            // add exclusion to parent carbon
            nbb.atom(start - 1).exclusion().insert(start + i);
            // set cg
            if (i == number - 1) {
              atn.setChargeGroup(cg);
            } else {
              atn.setChargeGroup(0);
            }
            ostringstream n;
            n << "H" << unique_id;
            if (number > 1) n  << (i+1);
            atn.setName(n.str());
            nbb.addAtom(atn);
          }

          // Inherit H exclusions from carbon
          AtomTopology at(nbb.atom(start - 1));
          Exclusion e;
          for (int i = 0; i < at.exclusion().size(); i++) {
            e.insert(at.exclusion().atom(i));
          }
          for (int i = 0; i < number; ++i) {
            e.erase(start + i);
            nbb.atom(start + i).setExclusion(e);
          }

          //now, for all the following atoms, we have to adapt the exclusions
          //without checking
          for (int i = start; i < obb.numAtoms(); i++) {
            AtomTopology at(obb.atom(i));
            Exclusion e;
            for (int j = 0; j < at.exclusion().size(); j++)
              e.insert(at.exclusion().atom(j) + number);
            at.setExclusion(e);

            nbb.addAtom(at);
          }

          // Now, we do the bonds
          vector<int> nghs;
          {
            BondIterator bi(obb);
            for (; bi; ++bi) {
              int a[2];
              for (int i = 0; i < 2; i++) {
                a[i] = bi()[i];
                if (a[i] >= start)a[i] += number;
              }
              if (a[0] == start - 1) nghs.push_back(a[1]);
              else if (a[1] == start - 1) nghs.push_back(a[0]);
              Bond b(a[0], a[1]);
              b.setType(bi().type());
              nbb.addBond(b);
            }
            for (int i = 0; i < number; ++i) {
              Bond b(start - 1, start + i);
              b.setType(2);
              nbb.addBond(b);
            }
          }
          {
            AngleIterator bi(obb);
            for (; bi; ++bi) {
              int a[3];
              for (int i = 0; i < 3; i++) {
                a[i] = bi()[i];
                if (a[i] >= start)a[i] += number;
              }
              Angle b(a[0], a[1], a[2]);
              b.setType(bi().type());
              nbb.addAngle(b);
            }
            // Add angles
            for (vector<int>::const_iterator it = nghs.begin(); it != nghs.end(); ++it) {
              for (int i = 0; i < number; ++i) {
                Angle b(*it, start - 1, start + i);
                b.setType(9);
              nbb.addAngle(b);
              }
            }
            for (int i = 0; i < number; ++i) {
              for (int j = i + 1; j < number; ++j) {
                Angle b(start + i, start - 1, start + j);
                b.setType(9);
                nbb.addAngle(b);
              }
            }
          }
          {
            ImproperIterator bi(obb);
            for (; bi; ++bi) {
              int a[4];
              for (int i = 0; i < 4; i++) {
                a[i] = bi()[i];
                if (a[i] >= start)a[i] += number;
              }
              Improper b(a[0], a[1], a[2], a[3]);
              b.setType(bi().type());
              nbb.addImproper(b);
            }
          }
          {
            DihedralIterator bi(obb);
            for (; bi; ++bi) {
              int a[4];
              for (int i = 0; i < 4; i++) {
                a[i] = bi()[i];
                if (a[i] >= start)a[i] += number;
              }
              Dihedral b(a[0], a[1], a[2], a[3]);
              b.setType(bi().type());
              nbb.addDihedral(b);
            }
          }
          {
            LJExceptionIterator li(obb);
            for(; li; ++li) {
              int a[2];
              for(int i = 0; i < 2; ++i) {
                a[i] = li()[i];
                if (a[i] >= start) a[i] += number;
              }
              LJException l(a[0], a[1]);
              l.setType(li().type());
              l.indicate() = li().indicate();
              for(set<int>::const_iterator it = li().cond().begin(),
                      to = li().cond().end(); it != to; ++it) {
                if (*it > start)
                  l.addCond(*it + number);
                else
                  l.addCond(*it);
              }
              nbb.addLJException(l);
            }
          }
          obb = nbb;
        }
      }
    }
    //Now we need the output, but there is no outBuildingBlock yet!
    cout.precision(5);
    cout << "# This is a prepared building block" << endl;
    cout << "# It is based on " << obb.resName() << endl;
    cout << "# Every united carbon atom was replaced by bare" << endl;
    cout << "# carbon and corresponding number of hydrogen atoms" << endl;
    cout << "# Bonds and angles were added as well" << endl;

    OutBuildingBlock outbb(cout);
    obb.setResName(new_name);
    outbb.writeSingle(obb, OutBuildingBlock::BBTypeSolute);
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}


