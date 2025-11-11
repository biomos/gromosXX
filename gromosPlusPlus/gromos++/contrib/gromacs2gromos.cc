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
 * @file gromacs2gromos.cc
 * converts gromacs topologies and trajectories
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor gromacs2gromos
 * @section gromacs2gromos converts GROMACS topologies and trajectories
 * @author @ref ns
 * @date 06.07.2010
 *
 * Program gromacs2gromos can be used to convert a GROMACS topology to GROMOS.
 * The topology (GROMACS .trn file) is read from the file given by \@gmxtop,
 * converted to GROMOS and is written to \@outtop. Please note that not all
 * features of the topology can be converted to GROMOS due to different
 * functionality and philosophies of the GROMACS and GROMOS codes.
 *
 * Supported topological information:
 * - atoms (mass, name, charge, type, exclusions, 1,4 interactions)
 * - residue numbers and names
 * - bonds, angles, imporper and proper dihedrals
 * - lennard-jones parameters (including 1,4 interactions)
 * - atom type names
 * - physical constants
 * - constraints
 * - settle constraint
 * - solvent (must be called SOL)
 *
 * Explicitly NOT supported information:
 * - perturbation
 * - other force field terms (you will be warned about them)
 *
 * In addition this program can also convert GROMACS trajectories (*.trr files)
 * given by \@gmxtraj and writes them either to a file (\@outtraj) or the standard
 * output. The latter is suggested as it can be piped though gzip on fly. Only
 * the timestep, coordinate and box information are extracted.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@gmxtopo</td><td>&lt;GROMACS topology file (*.trn)&gt; </td></tr>
 * <tr><td> \@outtopo</td><td>&lt;the GROMOS topology written&gt;] </td></tr>
 * <tr><td>[\@gmxtraj</td><td>&lt;GROMACS trajectory file (*.trr)&gt;]</td></tr>
 * <tr><td>[\@outtraj</td><td>&lt;the GROMOS trajectory written&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 gromacs2gromos
     @gmxtopo  gromacs_input.trn
     @outtopo  gromos_topo.top
     @gmxtraj  gromacs_traj.trr    | gzip > gromos_traj.trc.gz
   @endverbatim

 * <hr>
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"


#include "../config.h"
#include "../src/gromos/Exception.h"

#ifdef HAVE_GMX
extern "C" {
#include <gromacs/typedefs.h>
#include <gromacs/tpxio.h>
#include <gromacs/mtop_util.h>
#include <gromacs/physics.h>
}
#include <gromacs/smalloc.h>
extern "C" {
#include <gromacs/trnio.h>
}
#endif

using namespace args;
using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;

int main(int argc, char** argv) {

#ifdef HAVE_GMX
  Argument_List knowns;
  knowns << "gmxtopo" << "outtopo" << "gmxtraj" << "outtraj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@gmxtopo     <gromacs checkpoint/topology file>\n";
  usage += "\t@outtopo     <the topology written>\n";
  usage += "\t[@gmxtraj    <gromacs trajectory file>]\n";
  usage += "\t[@outtraj    <the trajectory written, or to stdout>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    real t;
    t_state state;
    t_tpxheader tpx;
    t_inputrec ir;
    gmx_mtop_t mtop;
    t_topology top;
    char * fn = const_cast<char*>(args["gmxtopo"].c_str());

    read_tpxheader(fn, &tpx, TRUE, NULL, NULL);
    read_tpx_state(fn, tpx.bIr  ? &ir : NULL, &state,NULL, tpx.bTop ? &mtop: NULL);
    top = gmx_mtop_t_to_t_topology(&mtop);

    LinearTopology linTop;
    GromosForceField ggf;
    
    // add the atoms and their basic parameters
    cerr << "Adding atoms...";
    for(int i = 0; i < top.atoms.nr; ++i) {
      AtomTopology atom;
      ostringstream name;
      name << *(top.atoms.atomname[i]);
      atom.setName(name.str());
      atom.setCharge(top.atoms.atom[i].q);
      atom.setIac(top.atoms.atom[i].type);
      atom.setMass(top.atoms.atom[i].m);
      linTop.addAtom(atom);
      linTop.setResNum(i, top.atoms.resinfo[top.atoms.atom[i].resind].nr);
      //cerr << "residue: " << top.atoms.atom[i].resnr << endl;
    }
    // add the charge groups
    int start = 0, end = 0;
    for(int i = 0; i < top.cgs.nr; ++i) {
      end = top.cgs.index[i+1];
      //cerr << "cg["<< i << "] = {" << start << ".." << end-1 << "}" <<endl;
      for(int j = start; j < end; ++j) {
        if (j == end - 1) 
          linTop.atoms()[j].setChargeGroup(1);
        else
          linTop.atoms()[j].setChargeGroup(0);
      }
      start = end;
    }
    // add the exclusions
    start = 0; end = 0;
    for(int i = 0; i < top.excls.nr; ++i) {
      end = top.excls.index[i+1];
      //cerr << "cg["<< i << "][" << start << ".." << end-1 << "] = ";
      for(int j = start; j < end; ++j) {
        //cerr << top.excls.a[j] << " ";
        int ex = top.excls.a[j];
        if (i < ex)
          linTop.atoms()[i].exclusion().insert(ex);
        else if (i > ex)
          linTop.atoms()[ex].exclusion().insert(i);
      }
      //cerr << endl;
      start = end;
    }
    // add the 14 exclusions
    t_iatom *iatoms = top.idef.il[F_LJ14].iatoms;
    for(int i = 0; i < top.idef.il[F_LJ14].nr; i+=3) {
      ++iatoms; // don't need the type
      int a = *(iatoms++);
      int b = *(iatoms++);
      //cerr << "LJ14: " << a << " " << b << endl;
      if (a < b) {
        linTop.atoms()[a].exclusion14().insert(b);
        linTop.atoms()[a].exclusion().erase(b);
      } else if (a > b) {
        linTop.atoms()[b].exclusion14().insert(a);
        linTop.atoms()[b].exclusion().erase(a);
      }
    }
    cerr << " done." << endl;

    // add the residues
    cerr << "Adding " << top.atoms.nres << " residues...";
    for(int i = 0; i < top.atoms.nres; ++i) {
      ostringstream name;
      name << *(top.atoms.resinfo[top.atoms.atom[i].resind].name);
      linTop.setResName(i, name.str());
    }
    cerr << " done." << endl;

    // the bond types
    map<int, int> ffType;
    map<int, pair<int, int> > settleType;
    cerr << "Adding bonds...";
    int currentBondType = 0;
    for(int i = 0; i < top.idef.ntypes; ++i) {
      if (top.idef.functype[i] == F_BONDS) {
        BondType bt(currentBondType, 0.0, top.idef.iparams[i].harmonic.krA, top.idef.iparams[i].harmonic.rA);
        ffType[i] = currentBondType;
        ++currentBondType;
        ggf.addBondType(bt);
      } else if (top.idef.functype[i] == F_G96BONDS) {
        BondType bt(currentBondType, top.idef.iparams[i].harmonic.krA, 0.0, sqrt(top.idef.iparams[i].harmonic.rA));
        ffType[i] = currentBondType;
        ++currentBondType;
        ggf.addBondType(bt);
      } else if (top.idef.functype[i] == F_SETTLE) {
        //cerr << "settle types " << i << " " << top.idef.iparams[i].settle.doh << " " << top.idef.iparams[i].settle.dhh << endl;
        //cerr << "    will be known as: " << currentBondType << " and " << currentBondType + 1 << endl;
        BondType doh(currentBondType, 0.0, 0.0, top.idef.iparams[i].settle.doh);
        BondType dhh(currentBondType+1, 0.0, 0.0, top.idef.iparams[i].settle.dhh);
        settleType[i] = pair<int, int>(currentBondType, currentBondType+1);
        currentBondType += 2;
        ggf.addBondType(doh);
        ggf.addBondType(dhh);
      } else if (top.idef.functype[i] == F_CONSTR) {
        BondType bt(currentBondType, 0.0, 0.0, top.idef.iparams[i].constr.dA);
        ffType[i] = currentBondType;
        ++currentBondType;
        ggf.addBondType(bt);
      }
    }

    // the bonds
    iatoms = top.idef.il[F_BONDS].iatoms;
    for(int i = 0; i < top.idef.il[F_BONDS].nr; i+=3) {
      int t = *(iatoms++); 
      int a = *(iatoms++);
      int b = *(iatoms++);
      //cerr << "bond: " << t << " " << a << " " << b << endl;
      Bond bond(a,b);
      bond.setType(ffType[t]);
      linTop.addBond(bond);
    }
    iatoms = top.idef.il[F_G96BONDS].iatoms;
    for(int i = 0; i < top.idef.il[F_G96BONDS].nr; i+=3) {
      int t = *(iatoms++); 
      int a = *(iatoms++);
      int b = *(iatoms++);
      //cerr << "bond: " << t << " " << a << " " << b << endl;
      Bond bond(a,b);
      bond.setType(ffType[t]);
      linTop.addBond(bond);
    }

    iatoms = top.idef.il[F_SETTLE].iatoms;
    for(int i = 0; i < top.idef.il[F_SETTLE].nr; i+=2) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      //cerr << "settle: " << t << " " << a << endl;
      Bond oh1bond(a,a+1);
      oh1bond.setType(settleType[t].first);
      linTop.addBond(oh1bond);
      Bond oh2bond(a,a+2);
      oh2bond.setType(settleType[t].first);
      linTop.addBond(oh2bond);
      Bond hhbond(a+1,a+2);
      hhbond.setType(settleType[t].second);
      linTop.addBond(hhbond);
    }
    iatoms = top.idef.il[F_CONSTR].iatoms;
    for(int i = 0; i < top.idef.il[F_CONSTR].nr; i+=3) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      //cerr << "constraint: " << t << " " << a << " " << b << endl;
      Bond bond(a,b);
      bond.setType(ffType[t]);
      linTop.addBond(bond);
    }
    cerr << " done." << endl;

    cerr << "Adding angles...";
    int currentAngleType = 0;
    for(int i = 0; i < top.idef.ntypes; ++i) {
      if (top.idef.functype[i] == F_ANGLES) {
        AngleType bt(currentAngleType, 0.0, top.idef.iparams[i].harmonic.krA, top.idef.iparams[i].harmonic.rA);
        ffType[i] = currentAngleType;
        ++currentAngleType;
        ggf.addAngleType(bt);
      } else if (top.idef.functype[i] == F_G96ANGLES) {
        AngleType bt(currentAngleType, top.idef.iparams[i].harmonic.krA, 0.0, 180.0*acos(top.idef.iparams[i].harmonic.rA)/M_PI);
        ffType[i] = currentAngleType;
        ++currentAngleType;
        ggf.addAngleType(bt);
      }
    }

    // the angles
    iatoms = top.idef.il[F_ANGLES].iatoms;
    for(int i = 0; i < top.idef.il[F_ANGLES].nr; i+=4) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      int c = *(iatoms++);
      //cerr << "angle: " << t << " " << a << " " << b << " " << c << endl;
      Angle angle(a,b,c);
      angle.setType(ffType[t]);
      linTop.addAngle(angle);
    }
    iatoms = top.idef.il[F_G96ANGLES].iatoms;
    for(int i = 0; i < top.idef.il[F_G96ANGLES].nr; i+=4) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      int c = *(iatoms++);
      //cerr << "angle: " << t << " " << a << " " << b << " " << c << endl;
      Angle angle(a,b,c);
      angle.setType(ffType[t]);
      linTop.addAngle(angle);
    }
    cerr << " done." << endl;
    
    cerr << "Adding impropers...";
    int currentImproperType = 0;
    for(int i = 0; i < top.idef.ntypes; ++i) {
      if (top.idef.functype[i] == F_IDIHS) {
        ImproperType bt(currentImproperType, 
                top.idef.iparams[i].harmonic.krA * M_PI * M_PI / 180.0 / 180.0,
                top.idef.iparams[i].harmonic.rA);
        ffType[i] = currentImproperType;
        ++currentImproperType;
        ggf.addImproperType(bt);
      }
    }

    // the impropers
    iatoms = top.idef.il[F_IDIHS].iatoms;
    for(int i = 0; i < top.idef.il[F_IDIHS].nr; i+=5) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      int c = *(iatoms++);
      int d = *(iatoms++);
      //cerr << "improper: " << t << " " << a << " " << b << " " << c << " " << d << endl;
      Improper improper(a,b,c,d);
      improper.setType(ffType[t]);
      linTop.addImproper(improper);
    }
    cerr << " done." << endl;

    cerr << "Adding dihedrals...";
    int currentDihedralType = 0;
    for(int i = 0; i < top.idef.ntypes; ++i) {
      if (top.idef.functype[i] == F_PDIHS) {
        DihedralType bt(currentDihedralType, 
                        top.idef.iparams[i].pdihs.cpA, 
                        cos(M_PI * top.idef.iparams[i].pdihs.phiA / 180.0), 
                        top.idef.iparams[i].pdihs.phiA,
                        top.idef.iparams[i].pdihs.mult); 
        ffType[i] = currentDihedralType;
        ++currentDihedralType;
        ggf.addDihedralType(bt);
      }
    }

    // the dihedrals
    iatoms = top.idef.il[F_PDIHS].iatoms;
    for(int i = 0; i < top.idef.il[F_PDIHS].nr; i+=5) {
      int t = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      int c = *(iatoms++);
      int d = *(iatoms++);
      //cerr << "dihedral: " << t << " " << a << " " << b << " " << c << " " << d << endl;
      Dihedral dihedral(a,b,c,d);
      dihedral.setType(ffType[t]);
      linTop.addDihedral(dihedral);
    }
    cerr << " done." << endl;

    cerr << "Parsing system...";
    System sys(linTop.parse());
    cerr << " done." << endl;

    // set the physical constants for the units
    ggf.setBoltz(BOLTZ);
    ggf.setHbar(PLANCK / (2.0 * M_PI));
    ggf.setFpepsi(ONE_4PI_EPS0);
    ggf.setSpdl(SPEED_OF_LIGHT);
    ggf.setForceField("GMX");

    // the LJ matrix
    cerr << "Creating LJ matrix...";
    for(int i = 0, k = 0; i < top.idef.atnr; ++i) {
      for(int j = 0; j < top.idef.atnr; ++j, ++k) {
        AtomPair p(i,j);
        LJType lj(top.idef.iparams[k].lj.c12, top.idef.iparams[k].lj.c6);
        ggf.setLJType(p, lj);
      }
    }

    // the 1,4 entries
    iatoms = top.idef.il[F_LJ14].iatoms;
    for(int i = 0; i < top.idef.il[F_LJ14].nr; i+=3) {
      int k = *(iatoms++);
      int a = *(iatoms++);
      int b = *(iatoms++);
      AtomPair p(top.atoms.atom[a].type, top.atoms.atom[b].type);
      LJType lj(ggf.ljType(p));
      lj.cs12() = top.idef.iparams[k].lj14.c12A;
      lj.cs6() = top.idef.iparams[k].lj14.c6A;
      ggf.setLJType(p, lj);
    }

    // atom type names
    map<int,string> atomTypeNames;
    for(int i = 0; i < top.atoms.nr; ++i) {
      ostringstream name;
      name << *(top.atoms.atomtype[i]);
      atomTypeNames[top.atoms.atom[i].type] = name.str();
    }
    for(int i = 0; i < top.idef.atnr; ++i) {
      if (atomTypeNames.count(i))
        ggf.addAtomTypeName(atomTypeNames[i]);
      else
        ggf.addAtomTypeName("EMPTY");
    }
    cerr << " done." << endl;

    // detect solvent
    System outSys;
    int solMol = 0, solMolIndex = -1;
    int solAtoms = 0;
    for(int m = 0; m < sys.numMolecules(); ++m) {
      if (sys.mol(m).topology().resName(0) == "SOL") {
        ++solMol;
        solMolIndex = m;
        solAtoms += sys.mol(m).topology().numAtoms();
      } else {
        outSys.addMolecule(sys.mol(m));
      }
    }

    // issue some warnings about other FF term
    for(int i = 0; i < top.idef.ntypes; ++i) {
      switch(top.idef.functype[i]) {
        case F_BONDS :
        case F_G96BONDS :
        case F_ANGLES :
        case F_G96ANGLES :
        case F_PDIHS :
        case F_IDIHS :
        case F_LJ14 :
        case F_LJ :
        case F_CONSTR :
        case F_SETTLE :
          break;
        default:
          cerr << "WARNING: Dont know about interaction " << i + 1 << " of type "
                  << interaction_function[top.idef.functype[i]].name << " ("
                  << interaction_function[top.idef.functype[i]].longname << ")." << endl;
      }
    }

    if (solMol != 0) {
      cerr << "Adding " << solAtoms << " solvent atoms" << endl;
      SolventTopology soltop;
      for(int a = 0; a < sys.mol(solMolIndex).topology().numAtoms(); ++a) {
        soltop.addAtom(sys.mol(solMolIndex).topology().atom(a));
      }
      for(BondIterator b(sys.mol(solMolIndex).topology()); b; ++b) {
        Constraint c(b()[0], b()[1]);
        c.setDist(ggf.bondType(b().type()).b0());
        //cerr << "constraint type: " << b().type()<< endl;
        soltop.addConstraint(c);
      }
      Solvent sol(soltop);
      outSys.addSolvent(sol);
    } else {
      // add a dummy solvnet.
      cerr << "No solvent found. Adding dummy solvent." << endl;
      AtomTopology solatom;
      solatom.setCharge(0.0);
      solatom.setMass(1.008);
      solatom.setIac(0);
      solatom.setName("DUM");
      SolventTopology soltop;
      soltop.addAtom(solatom);
      outSys.addSolvent(Solvent(soltop));
    }

    sys = outSys;

    int offset = 0;
    for(int m = 0; m < sys.numMolecules(); ++m) {
      offset += sys.mol(m).topology().numAtoms();
      sys.addTemperatureGroup(offset);
      sys.addPressureGroup(offset);
    }

    cerr << "The topology is written to " << args["outtopo"] << "...";
    ofstream of(args["outtopo"].c_str());
    if (!of.is_open())
      throw gromos::Exception(argv[0], "Cannot write the output topology.");
    
    OutTopology out(of);
    ostringstream title;
    title << *(top.name);
    out.setTitle(title.str());

    // set the H-mass in order to sort the bonded interactions
    for(int m = 0; m < sys.numMolecules(); ++m)
      sys.mol(m).topology().setHmass(1.008);

    out.write(sys, ggf);
    of.close();
    cerr << " done." << endl;

    if (args.count("gmxtraj") == -1)
      return 1;

    ostream * outFile = &cout;
    if (args.count("outtraj") != -1) {
      of.open(args["outtraj"].c_str());
      if (!of.is_open())
        throw gromos::Exception(argv[0], "Cannot write the output trajectory.");
      outFile = &of;
    }
    OutG96S oc(*outFile);
    oc.select("ALL");
    oc.writeTitle("Converted GROMACS trajectory...");

    // convert the traj
    std::cerr << "Converting the trajectories...";
    for(Arguments::const_iterator it = args.lower_bound("gmxtraj"), to = args.upper_bound("gmxtraj");
            it != to; ++it) {      
      t_trnheader trn;
      gmx_bool bOK;
      t_fileio * fpread = open_trn(const_cast<char*>(it->second.c_str()), const_cast<char*>("r"));
      int nframe = 0;
      while (fread_trnheader(fpread, &trn, &bOK)) {
        rvec *x = new rvec[trn.natoms];
        matrix box;
        if (fread_htrn(fpread, &trn,
            trn.box_size ? box : NULL,
            trn.x_size ? x : NULL,
            NULL, NULL)) {
          oc.writeTimestep(trn.step, trn.t);
          int i = 0;
          for(int m = 0; m < sys.numMolecules(); ++m) {
            sys.mol(m).initPos();
            for(int a = 0; a < sys.mol(m).numAtoms(); ++a, ++i) {
              sys.mol(m).pos(a) = Vec(x[i][0], x[i][1], x[i][2]);
            }
          }
          for(int s = 0; s < solAtoms; ++s, ++i) {
            sys.sol(0).addPos(Vec(x[i][0], x[i][1], x[i][2]));
          }

          sys.box() = Box(
                  Vec(box[0][0],box[1][0],box[2][0]),
                  Vec(box[0][1],box[1][1],box[2][1]),
                  Vec(box[0][2],box[1][2],box[2][2]));
          sys.hasBox = true;
          sys.box().setNtb(Box::triclinic);
            
          oc << sys;
        }
        
        delete[] x;
        nframe++;
      }
    } // for traj files
    cerr << " done." << endl;
    
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
#else
  cerr << "This program can only be used if GROMOS++ is linked with the"
          " GROMACS libraries" << endl << endl
          << "Recompile with 'configure --with-gromacs=/path/to/gromacs'." << endl;
  return 1;
#endif
}

