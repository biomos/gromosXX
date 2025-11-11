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
 * @file add_atom.cc
 * Insert atoms in a building block
 */

/**
 * @page contrib Contrib Documentation
 *
 * @anchor add_atom
 * @section add_atom Insert atoms in a building block
 * @author @ref co
 * @date 7. 6. 2007
 *
 * This program allows the user to modify existing building blocks by 
 * inserting additional atoms in the building block. The user specifies after 
 * which atom the insertion is to take place (\@start) and how many atoms are 
 * to be inserted (\@number).
 *
 * The program does not add any parameters for the new atoms, nor does it 
 * connect them through any covalent interactions. It only adopts all 
 * references to atoms after the insertion in the exclusions, bond, bond
 * angles, impropers and dihedrals, thereby preparing the building block for
 * a manual adaptation.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@file</td><td>&lt;molecular building block file&gt; </td></tr>
 * <tr><td> \@build</td><td>&lt;name of the building block to modify&gt; </td></tr>
 * <tr><td> \@start</td><td>&lt;index number of the atom after which the insertion should take place&gt; </td></tr>
 * <tr><td> \@number</td><td>&lt;number of atoms to insert&gt; </td></tr>
 * </table>
 *
 *
 * Example 1:
 * @verbatim
  add_atom
    @file    53a6.mtb
    @build   UREA
    @start   4
    @number  2
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <set>
#include <string>

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
  knowns << "file" << "build" << "start" << "number";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@file <mtb-file>\n";
  usage += "\t@build <buildingblock name>\n";
  usage += "\t@start <after which atom do you want to add>\n";
  usage += "\t@number <number of atoms to add>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    InBuildingBlock ibb(args["file"]);
    BuildingBlock mtb(ibb.building());

    int index = mtb.findBb(args["build"]);
    BbSolute obb(mtb.bb(index - 1));
    BbSolute nbb;

    int start = args.getValue<int>("start");
    int number = args.getValue<int>("number");

    if (start < 0 || number < 0)
      throw gromos::Exception(argv[0], "Start or number are negativ.");

    //prepare an empty atom
    AtomTopology atn;
    atn.setIac(0);
    atn.setMass(0);
    atn.setCharge(0.0);
    atn.setChargeGroup(0);
    atn.setName("New");

    //copy obb to nbb, but stop at start
    //first the previous exclusions
    for (int i = 0; i < obb.numPexcl(); i++) {
      Exclusion e;
      for (int j = 0; j < obb.pexcl(i).size(); j++)
        if (obb.pexcl(i).atom(j) < start)
          e.insert(obb.pexcl(i).atom(j));
        else
          e.insert(obb.pexcl(i).atom(j) + number);

      nbb.addPexcl(e);
    }

    //and the name
    nbb.setResName(obb.resName());

    //now, the atoms bis und mit start, pay attention to the exclusions.
    for (int i = 0; i < start; i++) {
      AtomTopology at(obb.atom(i));
      Exclusion e;
      for (int j = 0; j < at.exclusion().size(); j++)
        if (at.exclusion().atom(j) < start)
          e.insert(at.exclusion().atom(j));
        else
          e.insert(at.exclusion().atom(j) + number);
      at.setExclusion(e);
      nbb.addAtom(at);
    }

    //add number atoms
    for (int i = 0; i < number; i++)
      nbb.addAtom(atn);
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
    {
      BondIterator bi(obb);
      for (; bi; ++bi) {
        int a[2];
        for (int i = 0; i < 2; i++) {
          a[i] = bi()[i];
          if (a[i] >= start)a[i] += number;
        }
        Bond b(a[0], a[1]);
        b.setType(bi().type());
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
    //Now we need the output, but there is no outBuildingBlock yet!
    cout.precision(5);
    cout << "# This is a prepared building block" << endl;
    cout << "# It is based on " << nbb.resName() << endl;
    cout << "# I have added " << number << " atoms, starting after atom "
            << start << endl;
    cout << "# I only corrected the exclusions, bonds, angles etc." << endl;
    cout << "# of the EXISTING atoms. You still have to do the following: "
            << endl;
    cout << "# 1. Adapt the exclusions of atoms before atom " << start + 1
            << " to exclude the new atoms" << endl;
    cout << "# 2. Adapt all parameters of the new atoms, including exclusions"
            << endl;
    cout << "# 3. Add all bonds, angles, etc. involving any new atoms"
            << endl;
    cout << "# 4. Possibly adapt bonds, angles, etc around the new atoms"
            << endl;
    cout << "#" << endl;

    OutBuildingBlock outbb(cout);
    outbb.writeSingle(nbb, OutBuildingBlock::BBTypeSolute);
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}


