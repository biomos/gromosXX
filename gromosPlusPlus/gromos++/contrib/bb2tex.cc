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
 * @file bb2tex.cc
 * creates latex tables from a building block
 */
/**
 * @page contrib Contrib program documentation
 *
 * @anchor bb2tex
 * @section bb2tex convert a building block to latex
 * @author @ref ns
 * @date 01-07-2009
 *
 * The program bb2tex takes a building block and writes tables for the atoms,
 * bonds, dihedrals and improper dihedrals.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@build</td><td>&lt; molecular building block file &gt; </td></tr>
 * <tr><td> \@seq</td><td>&lt; block of interest &gt; </td></tr>
 */


#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;

int main(int argc, char *argv[]) {

  Argument_List knowns;
  knowns << "build" << "seq";

  string usage = string("# ") + argv[0];
  usage += "\n\t@build <mtb-file>\n";
  usage += "\t@seq <building block>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read in the building block file
    BuildingBlock mtb;
    {
      Arguments::const_iterator iter = args.lower_bound("build"),
              to = args.upper_bound("build");
      for (; iter != to; ++iter) {
        InBuildingBlock ibb(iter->second);
        mtb.addBuildingBlock(ibb.building());
      }
    }

    // loop over blocks.
    for(Arguments::const_iterator iter = args.lower_bound("seq"),
            to = args.upper_bound("seq"); iter != to; ++iter) {
      const string & name = iter->second;
      int index = mtb.findBb(name);
      if (index == 0) {
        ostringstream msg;
        msg << "Building block " << name << " was not found.";
        throw gromos::Exception(argv[0], msg.str());
      }

      BbSolute const * p_bb = NULL;
      int last_few=0;
      if (index < 0) {
        cerr << name << " is an end group." << endl;
        p_bb = &mtb.be(abs(index)-1);
        last_few=p_bb->rep();
      } else {
        p_bb = &mtb.bb(index-1);
        last_few=p_bb->numPexcl();
      }

      const BbSolute & block = *p_bb;
      ofstream atoms_file((name + "_atoms.tex").c_str());
      if (!atoms_file.is_open())
        throw gromos::Exception(argv[0], "Cannot write files.");

      atoms_file << "\\gromostabular{" << endl
              << "Seq. & Name & IAC & Mass & Charge & Exclusions" << endl
              << "}" << endl;
      // preceding exlusions
      for(int i=0; i < block.numPexcl(); i++){
        int index = i - block.numPexcl() + 1;
	atoms_file << setw(5) << index << "&      &       &       &             &";
	for(int j=0; j< block.pexcl(i).size();j++)
	  atoms_file << " " << block.pexcl(i).atom(j)+1;
	atoms_file << " \\n" << endl;
      }
      if (block.numPexcl())
        atoms_file << "\\hline" << endl;

      // atoms
      for(int i = 0; i < block.numAtoms(); ++i) {
        atoms_file << setw(5) << i+1 << " & ";

        atoms_file.setf(ios::left, ios::adjustfield);

        atoms_file << setw(4) << block.atom(i).name() << " & ";
        atoms_file.setf(ios::fixed, ios::adjustfield);
        atoms_file.precision(5);
        atoms_file.setf(ios::fixed, ios::floatfield);

        atoms_file << setw(5) << block.atom(i).iac() + 1 << " & "
                << setw(5) << int(block.atom(i).mass()) + 1 << " & "
                << setw(11) << block.atom(i).charge() << " &";

        if (i < block.numAtoms() - last_few) {
          for (int j = 0; j < block.atom(i).exclusion().size(); j++) {
            atoms_file << " " << block.atom(i).exclusion().atom(j) + 1;
          }
        }
        if (i != block.numAtoms() - 1) {
          if (block.atom(i).chargeGroup() == 1)
            atoms_file << " \\n" << endl << "\\hline" << endl;
          else
            atoms_file << " \\n" << endl;
        }
      }
      atoms_file.close();

      // bonds
      ofstream bonds_file((name + "_bonds.tex").c_str());
      if (!bonds_file.is_open())
        throw gromos::Exception(argv[0], "Cannot write files.");

      bonds_file << "\\gromostabular{" << endl
              << " I & J & Type" << endl
              << "}" << endl;

      for (BondIterator bi(block); bi; ++bi) {
        bonds_file << setw(5) << bi()[0] + 1 << " & "
                << setw(5) << bi()[1] + 1 << " & "
                << setw(5) << bi().type() + 1;
        if (!bi.last())
          bonds_file << " \\n" << endl;
      }
      bonds_file.close();

      // angles
      ofstream angles_file((name + "_angles.tex").c_str());
      if (!angles_file.is_open())
        throw gromos::Exception(argv[0], "Cannot write files.");

      angles_file << "\\gromostabular{" << endl
              << " I & J & K & Type" << endl
              << "}" << endl;

      for (AngleIterator ai(block); ai; ++ai) {
        angles_file << setw(5) << ai()[0] + 1 << " & "
                << setw(5) << ai()[1] + 1 << " & "
                << setw(5) << ai()[2] + 1 << " & "
                << setw(5) << ai().type() + 1;
        if (!ai.last())
          angles_file << " \\n" << endl;
      }
      angles_file.close();

      // dihedrals
      ofstream dihedrals_file((name + "_dihedrals.tex").c_str());
      if (!dihedrals_file.is_open())
        throw gromos::Exception(argv[0], "Cannot write files.");

      dihedrals_file << "\\gromostabular{" << endl
              << " I & J & K & L & Type" << endl
              << "}" << endl;

      for (DihedralIterator di(block); di; ++di) {
        dihedrals_file << setw(5) << di()[0] + 1 << " & "
                << setw(5) << di()[1] + 1 << " & "
                << setw(5) << di()[2] + 1 << " & "
                << setw(5) << di()[3] + 1 << " & "
                << setw(5) << di().type() + 1;
        if (!di.last())
          dihedrals_file << " \\n" << endl;
      }
      dihedrals_file.close();
      // impropers
      ofstream impropers_file((name + "_impropers.tex").c_str());
      if (!impropers_file.is_open())
        throw gromos::Exception(argv[0], "Cannot write files.");

      impropers_file << "\\gromostabular{" << endl
              << " I & J & K & L & Type" << endl
              << "}" << endl;

      for (ImproperIterator ii(block); ii; ++ii) {
        impropers_file << setw(5) << ii()[0] + 1 << " & "
                << setw(5) << ii()[1] + 1 << " & "
                << setw(5) << ii()[2] + 1 << " & "
                << setw(5) << ii()[3] + 1 << " & "
                << setw(5) << ii().type() + 1;
        if (!ii.last())
          impropers_file << " \\n" << endl;
      }
      impropers_file.close();
    }
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
