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
 * @file gch.cc
 * Generate coordinates for explicit hydrogens
 */

/**
 * @page programs Program Documentation
 *
 * @anchor gch
 * @section gch Generate coordinates for explicit hydrogens
 * @author @ref co
 * @date 11-6-07
 *
 * In the standard GROMOS force fields, part of the hydrogen atoms (polar,
 * aromatic) are explicitly treated, whereas other hydrogen atoms (aliphatic,
 * some aromatic) are implicitly treated by incorporation into the (carbon)atom
 * to which they are attached. Depending on the presence or absence of hydrogen
 * atom coordinates in a molecular configuration file, hydrogen atom
 * coordinates may have to be recalculated.
 *
 * Program gch calculates optimal positions for hydrogen atoms for which the 
 * connecting bond shows a relative deviation from the zero-energy value larger
 * than a user specified threshold. Coordinates for all hydrogen atoms that are
 * explicitly listed in the topology should already be contained in the
 * coordinate file. Program @ref pdb2g96 e.g. will include atoms for which no
 * coordinates were present in the pdb-file with coordinates set to zero. If
 * defined, gch uses topological information on bonds, bond angles and dihedral
 * angles to place hydrogen atoms at the optimal location. In cases where the
 * necessary angular parameters are not provided in the topology, gch uses
 * 109.5 degrees for tetrahedral centers and 120 degrees for planar centers.
 *
 * Eight types of geometries can be handled when generating hydrogen atom 
 * coordinates:
 * <ol>
 * <li> An atom (a) is bonded to one hydrogen (H) and one other heavy atom 
 *      (nh). A fourth atom (f) is searched for which is bounded to nh and 
 *      preferably is used to define the dihedral around the nh-a bond. The 
 *      coordinates of H are generated in such a way that the dihedral 
 *      f-nh-a-H is trans and that the angle nh-a-H and bond length a-H
 *      correspond to their minimum energy values.
 * <li> An atom (a) is bonded to one hydrogen (H) and two other heavy atoms
 *      (nh1 and nh2). The coordinates of H are generated to be in the plane 
 *      through nh1, nh2 and a, on the line bisecting the nh1-a-nh2 angle and
 *      with an a-H bond length corresponding to the minimum energy value in
 *      the topology, such that the nh1-a-H and nh2-a-H angles are larger than
 *      90 degrees.
 * <li> An atom (a) is bonded to two hydrogens (H1 and H2) and one other heavy
 *      atom (nh). A fourth atom (f) is searched for which is bounded to nh
 *      and preferably is used to define the dihedral around the nh-a bond. The
 *      coordinates of H1 are generated in such a way that the dihedral 
 *      f-nh-a-H1 is trans and that the angle nh-a-H1 and bond length a-H1 
 *      correspond to their minimum energy values. The coordinates of H2 are 
 *      generated to have the angles nh-a-H2 and H1-a-H2 as well as the bond
 *      length a-H2 at their minimum energy values. If this does not result in
 *      a planar configuration around a, the improper dihedral a-nh-H1-H2 will
 *      be positive.
 * <li> An atom (a) is bonded to three hydrogens (H1, H2 and H3) and one other
 *      heavy atom (nh). A fourth atom (f) is searched for wich is bounded to
 *      nh and preferably is used to define the dihedral around the nh-a bond.
 *      The coordinates of H1 are generated in such a way that the dihedral
 *      f-nh-a-H1 is trans and that the angle nh-a-H1 and bond length a-H1
 *      correspond to their minimum energy values. The coordinates of H2 are 
 *      such that the angles nh-a-H2 and H1-a-H2 and the bond length a-H2 are 
 *      at their minimum energy values, and the improper dihedral a-nh-H1-H2 is
 *      positive. The coordinates of H3 are such that the angles nh-a-H3 and
 *      H1-a-H3 and the bond length a-H3 are at their minimum energy values and
 *      the improper dihedral a-nh-H1-H3 has a negative value.
 * <li> An atom (a) is bonded to one hydrogen atom (H) and three other heavy
 *      atoms (nh1, nh2, nh3). The coordinates of H are generated along the line
 *      going through atom a and a point corresponding to the average position
 *      of nh1, nh2 and nh3, such that the bond length a-H is at its minimum
 *      energy value and the angles nh1-a-H, nh2-a-H and nh3-a-H are larger than
 *      90 degree.
 * <li> An atom (a) is bonded to two hydrogen atoms (H1 and H2) and two other
 *      heavy atoms (nh1 and nh2). The coordinates of H1and H2 are placed above
 *      and below the plane going through atoms nh1, nh2 and a, in such a way
 *      that the a-H1 and a-H2 bond lengths and the H1-a-H2 bond angle are at
 *      their minimum energy values. The improper dihedral angle a,nh1,nh2,H1
 *      will be positive.
 * <li> An atom (a) is bonded to two hydrogen atoms (H1 and H2), but to no
 *      heavy atoms. This is likely to be a (crystallographic) water molecule. 
 *      First a molecule is generated having the a-H1 aligned in the 
 *      z-direction and the a-H2 in the z-y plane with the angle H1-a-H2 and
 *      bond lengths a-H1 and a-H2 according to their minimum energy values. 
 *      This molecule is then rotated around x, y and z by three random angles.
 * <li> An atom (a) is bonded to four hydrogen atoms (H1, H2, H3 and H4), but
 *      to no heavy atoms. A molecule is generated with all bond lengths at 
 *      their minimum energy value, the a-H1 aligned in the z-directionm H2 in
 *      the x-z plane and H3 such that the improper a-H1-H2-H3 is positive and 
 *      H4 such that the improper a-H1-H2-H4 is negative. The complete molecule
 *      is then rotated by three random angles around x, y and z.
 * </ol>
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> \@tol</td><td>&lt;tolerance (default 0.1 %)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  gch
    @topo  ex.top
    @pos   exref.coo
    @tol   0.1
    @pbc   r
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Constraint.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Gch.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace utils;
using namespace gmath;
using namespace bound;



bool isnan(Vec v) {
  return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]);
}

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pos" << "tol" << "pbc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pos   <input coordinate file>\n";
  usage += "\t@tol   <tolerance (default 0.1 %)>\n";
  usage += "\t@pbc   <boundary type> <gather method>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology make system and force field
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // read in coordinates
    InG96 ic(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();
    System outSys(sys);

    // read in the accuracy
    double eps = args.getValue<double>("tol", false, 0.1) / 100.0;

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // gather the system!
    (*pbc.*gathmethod)();

    // a bit of an ugly hack for solvent molecules. The whole program has
    // been written for solutes, but sometimes we would have crystallographic
    // waters for which we want to do the same thing.
    if (sys.sol(0).numPos()) {
      MoleculeTopology mt;
      for (int a = 0; a < sys.sol(0).topology().numAtoms(); a++) {
        mt.addAtom(sys.sol(0).topology().atom(a));
        mt.setResNum(a, 0);
      }
      mt.setResName(0, "SOLV");

      ConstraintIterator ci(sys.sol(0).topology());
      for (; ci; ++ci) {
        gff.addBondType(BondType(gff.numBondTypes(), 1, ci().dist()));
        Bond b(ci()[0], ci()[1], false);
        b.setType(gff.numBondTypes() - 1);
        mt.addBond(b);
      }
      
      // add every solvent molecule as a solute
      int numSolvent = sys.sol(0).numPos() / sys.sol(0).topology().numAtoms();
      for (int i = 0; i < numSolvent; i++) {
        Molecule m(mt);
        m.initPos();
        for (int a = 0; a < mt.numAtoms(); a++) {
          m.pos(a) = sys.sol(0).pos(i * sys.sol(0).topology().numAtoms() + a);
        }
        sys.addMolecule(m);
      }
      // and remove the original solvent
      sys.sol(0).setNumPos(0);
    }

    // initialize two counters
    int replaced = 0, kept = 0;

    // loop over all atoms
    for (int m = 0; m < sys.numMolecules(); m++) {

      // flag the atoms with mass 1.008 as hydrogens
      sys.mol(m).topology().setHmass(1.008);
      for (int a = 0; a < sys.mol(m).numAtoms(); a++) {

        if (!sys.mol(m).topology().atom(a).isH()) {

        // divide into hydrogens and non-hydrogens
        vector<int> h;
        vector<int> nh;
        get_h_nh_neighbours(sys, gff, m, a, h, nh);

        // only continue if we have hydrogens
        int numH = h.size();
        int numNH = nh.size();
        int geom = get_geometry(numH, numNH);
        if (numH && !geom) {
          ostringstream os;
          os << "Unexpected geometry for hydrogen bound to atom: "
                  << m + 1 << ":" << a + 1 << endl;
          throw (gromos::Exception("gch", os.str()));
        }
        // we have to have a geometry (this means that there are hydrogens)
        // and a should not be a hydrogen itself. (in the case of H2O we have
        // e.g. H-H bonds. These are only treated via the oxygen.
        if (geom) {
          int r = generate_hcoordinates(sys, gff, m, a, h, nh, geom, eps);
          replaced += r;
          kept += (numH - r);
        }
        }
      }
    }

    bool found_nan=false;
    // fix the solvent hack
    int solventIndex = 0;
    for (int m = 0; m < sys.numMolecules(); ++m) {
      if (m < outSys.numMolecules()) {
        // solute
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a) {
          if (isnan(sys.mol(m).pos(a)))  {
            std::cerr << "Warning: "<<sys.mol(m).topology().resNum(a)+1 << " "
            << sys.mol(m).topology().resName(sys.mol(m).topology().resNum(a))             
            << " "  << sys.mol(m).topology().atom(a).name() 
            << " " <<v2s(sys.mol(m).pos(a)) << std::endl;
            found_nan=true;
            }
          outSys.mol(m).pos(a) = sys.mol(m).pos(a);
        }
      } else {
        // solvent
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++solventIndex) {
          if (isnan(sys.mol(m).pos(a))) {
            std::cerr << "Warning: "<<sys.mol(m).topology().resNum(a)+1 << " "
            << sys.mol(m).topology().resName(sys.mol(m).topology().resNum(a))             
            << " "  << sys.mol(m).topology().atom(a).name() 
            << " " <<v2s(sys.mol(m).pos(a)) << std::endl;
             found_nan=true;
           }
          outSys.sol(0).pos(solventIndex) = sys.mol(m).pos(a);
        }
      }
    }

    if (found_nan) {
      std::cerr << "WARNING: Some positions could not be generated (nan).\n  This can e.g. be due to missing or overlapping heteroatom positions.\n";   
    }
    
    ostringstream os;
    os <<"gch found " << replaced + kept << " hydrogen atoms in "
          << args["pos"] << endl;
    os << kept << " were within " << eps * 100 << "% of minimum energy bond "
          << "length" << endl;
    os << replaced << " were assigned new coordinates based on geometry";

    // now define an output stream and write the coordinates
    OutG96S oc;
    try{
      args.check("out", 1);
      ofstream fout(args["out"].c_str());
      oc.open(fout);
    }
    catch(const gromos::Exception &e){
      oc.open(cout);
    }
    oc.select("ALL");
    oc.writeTitle(os.str());
    oc << outSys;
    oc.close();

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
