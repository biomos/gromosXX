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

// Some functions needed by the program make_top
// several might be useful for later programs as well

#ifndef INCLUDED_MAKE_TOP
#define INCLUDED_MAKE_TOP

#include <iostream>
#include <string>
#include <sstream>

#include "../gromos/Exception.h"
#include "../gcore/BbSolute.h"
#include "../gcore/LinearTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Improper.h"
#include "../gcore/Dihedral.h"

using namespace std;
using namespace gcore;


void addSolute(gcore::LinearTopology &lt,
        BbSolute bb, int resnum, std::string resname, int rep, int nn);
int addBegin(gcore::LinearTopology &lt,
        BbSolute bb, int resnum);
void addEnd(gcore::LinearTopology &lt,
        BbSolute bb, int resnum);
void addCovEnd(gcore::LinearTopology &lt,
        BbSolute bb, int offset);

// as a little extra, a function to generate all 1,4-interactions based on the
// bonds and an ugly hack to put in Cysteine bonds
void setCysteines(gcore::LinearTopology &lt,
        int a, int b);
void setHeme(gcore::LinearTopology &lt,
        int a1, int a2, int b);
void prepareCyclization(gcore::LinearTopology &lt);
void cyclize(gcore::LinearTopology &lt);

// and a function that returns all atoms that are bonded to a set of atoms
// but have a lower number than offset
// This function is needed to determine the candidates for a bond
std::set<int> bondedAtoms(std::set<gcore::Bond> &bonds,
        std::set<int> atoms,
        int offset);

// ==========================================================================
//
// implementation:
//
// ===========================================================================

void addSolute(gcore::LinearTopology &lt,
        BbSolute bb, int resnum, std::string resname, int rep, int nn) {
    int na = lt.atoms().size();
    int strt = na - rep;
    int beg = 0;

    if (strt < 0) beg = -strt;

    for (int i = beg; i < rep; i++) {
        Exclusion e;
        for (int j = 0; j < bb.atom(i).exclusion().size(); j++)
            e.insert(bb.atom(i).exclusion().atom(j) + strt);
        lt.atoms()[strt + i].setExclusion(e);
    }

    //or, if rep=0, but we have preceding exclusions
    if (rep == 0) {
        int pexl = strt - bb.numPexcl();
        if (pexl < 0) throw gromos::Exception("addSolute",
                "Preceding exclusions, but no preceding atoms\n");

        /*
         *  Checking if there is a bond with a atom greater than the lt.numAtoms()
         * and bb.numPexcl() == 0;
         *
         */
        int offset = 1;
        int biggestBondedAtom = 0;
        for (std::set<Bond>::iterator it = lt.bonds().begin(); it != lt.bonds().end(); ++it) {
            Bond b(*it);
            if(b[0]+offset > biggestBondedAtom)
                biggestBondedAtom = b[0]+offset;
            if(b[1]+offset > biggestBondedAtom)
                biggestBondedAtom = b[1]+offset;
        }

        if(biggestBondedAtom > na && bb.numPexcl() == 0) {
            stringstream ss;
            ss << "Succeeding bonds at residue " << lt.resNames().back() << ", but no preceding exclusion on residue " << resname ;
            throw gromos::Exception("addSolute", ss.str());
        }


        for (int i = 0; i < bb.numPexcl(); i++) {
            Exclusion e;
            for (int j = 0; j < bb.pexcl(i).size(); j++)
                e.insert(bb.pexcl(i).atom(j) + strt);
            lt.atoms()[pexl + i].setExclusion(e);
        }
    }

    //finally, we add the atoms, but leave out the first rep
    for (int i = rep; i < bb.numAtoms(); i++) {
        Exclusion e;
        for (int j = 0; j < bb.atom(i).exclusion().size(); j++)
            e.insert(bb.atom(i).exclusion().atom(j) + strt);

        lt.addAtom(bb.atom(i));
        lt.atoms()[strt + i].setExclusion(e);
        lt.setResNum(strt + i, resnum);
    }
    lt.setResName(resnum, resname);

    // CG information
    for (int i = 0; i < bb.numAtoms(); i++) {
      lt.atoms()[strt+i].setCoarseGrained(bb.atom(i).isCoarseGrained());
      lt.atoms()[strt+i].setCGFactor(bb.atom(i).cg_factor());
    }

    // Polarisation
    for (int i = 0; i < bb.numAtoms(); i++) {
      lt.atoms()[strt+i].setPolarisable(bb.atom(i).isPolarisable());
      lt.atoms()[strt+i].setPolarisability(bb.atom(i).polarisability());
      lt.atoms()[strt+i].setCosCharge(bb.atom(i).cosCharge());
      lt.atoms()[strt+i].setDampingLevel(bb.atom(i).dampingLevel());
      lt.atoms()[strt+i].setDampingPower(bb.atom(i).dampingPower());
      lt.atoms()[strt+i].setPoloffsiteGamma(bb.atom(i).poloffsiteGamma());
      lt.atoms()[strt+i].setPoloffsiteI(bb.atom(i).poloffsiteI()+strt);
      lt.atoms()[strt+i].setPoloffsiteJ(bb.atom(i).poloffsiteJ()+strt);
    }

    //now, bonded interactions
    int offset = strt;

    //bonds
    BondIterator bi(bb);
    for (; bi; ++bi) {
        Bond b(bi()[0] + offset, bi()[1] + offset);
        b.setType(bi().type());

        //check if it exists already
        int found = 0;
        if (rep > 0)
            for (std::set<gcore::Bond>::const_iterator iter = lt.bonds().begin();
                    iter != lt.bonds().end(); ++iter)
                if ((*iter)[0] == b[0] && (*iter)[1] == b[1]) found = 1;

        //if(!found && b[0]>=nn && b[1]>=nn)
        // this leads to problems if we want a CYS1 / HIS1 as the first residue,
        // where we want a bond to a negative number until we parsed it... What
        // goes wrong if I leave the check out?
        if (!found)
            lt.addBond(b);
    }

    //dipole bonds
    BondDipoleIterator bdi(bb);
    for (; bdi; ++bdi) {
        Bond b(bdi()[0] + offset, bdi()[1] + offset);
        b.setType(bdi().type());

        //check if it exists already
        int found = 0;
        if (rep > 0)
            for (std::set<gcore::Bond>::const_iterator iter = lt.dipoleBonds().begin();
                    iter != lt.dipoleBonds().end(); ++iter)
                if ((*iter)[0] == b[0] && (*iter)[1] == b[1]) found = 1;

        //if(!found && b[0]>=nn && b[1]>=nn)
        // this leads to problems if we want a CYS1 / HIS1 as the first residue,
        // where we want a bond to a negative number until we parsed it... What
        // goes wrong if I leave the check out?
        if (!found)
            lt.addDipoleBond(b);
    }

    //angles
    AngleIterator ai(bb);
    for (; ai; ++ai) {
        Angle b(ai()[0] + offset, ai()[1] + offset, ai()[2] + offset);
        b.setType(ai().type());

        //check if it exists already
        int found = 0;
        if (rep > 0)
            for (std::set<gcore::Angle>::const_iterator iter = lt.angles().begin();
                    iter != lt.angles().end(); ++iter)
                if ((*iter)[0] == b[0] && (*iter)[1] == b[1] &&
                        (*iter)[2] == b[2]) found = 1;

        if (b[0] >= na || b[1] >= na || b[2] >= na)
            // if(!found&& b[0] >= nn && b[1] >= nn && b[2] >= nn)
            // what goes wrong if we leave out the check for negative numbers?
            if (!found)
                lt.addAngle(b);
    }

    //impropers
    ImproperIterator ii(bb);
    for (; ii; ++ii) {
        Improper b(ii()[0] + offset, ii()[1] + offset, ii()[2] + offset, ii()[3] + offset);
        b.setType(ii().type());

        //check if it exists already
        int found = 0;
        if (rep > 0)
            for (std::set<gcore::Improper>::const_iterator
                iter = lt.impropers().begin(); iter != lt.impropers().end(); ++iter)
                if ((*iter)[0] == b[0] && (*iter)[1] == b[1] &&
                        (*iter)[2] == b[2] && (*iter)[3] == b[3]) found = 1;

        if (b[0] >= na || b[1] >= na || b[2] >= na || b[3] >= na)
            // if(!found && b[0] >= nn && b[1] >= nn && b[2] >= nn && b[3] >= nn)
            // what goes wrong if we leave out the check for negative numbers?
            if (!found)
                lt.addImproper(b);
    }

    //dihedrals
    DihedralIterator di(bb);
    int counter = 0;

    // we do this in a seperate piece of code depending on rep
    if (rep == 0) {
        for (; di; ++di) {
            counter++;

            // find what is position -2 in the first dihedral
            int corr0 = offset;
            int corr1 = offset;
            int corr2 = offset;
            int corr3 = offset;
            // check for some heme properties first
            if (di()[3] == -2) {
                corr3 = 0;
            }

            if (di()[0] == -3) {
                corr0 = 0;
                for (std::set<gcore::Bond>::const_iterator iter = lt.bonds().begin();
                        iter != lt.bonds().end(); ++iter)
                    if ((*iter)[1] == di()[1] + offset) corr0 = (*iter)[0] + 3;
            }

            Dihedral b(di()[0] + corr0, di()[1] + corr1, di()[2] + corr2, di()[3] + corr3);
            b.setType(di().type());
            lt.addDihedral(b);
        }
    } else if (rep > 0) {
        // bad luck. We have to copy everything over several times
        std::vector<gcore::Dihedral> oldDihedral;
        for (std::set<gcore::Dihedral>::const_iterator iter = lt.dihedrals().begin();
                iter != lt.dihedrals().end(); ++iter)
            oldDihedral.push_back(*iter);
        lt.dihedrals().clear();

        // we start the same way as before
        for (; di; ++di) {
            counter++;

            // find what is position -2 in the first dihedral
            int corr = offset;
            if (di()[0] == -3) {
                for (std::set<gcore::Bond>::const_iterator iter = lt.bonds().begin();
                        iter != lt.bonds().end(); ++iter)
                    if ((*iter)[1] == di()[1] + offset) corr = (*iter)[0] + 3;
            }
            // or in this case di() = -1 might also be ill defined
            if (di()[0] == -2) {
                for (std::set<gcore::Bond>::const_iterator iter = lt.bonds().begin();
                        iter != lt.bonds().end(); ++iter)
                    if ((*iter)[1] == di()[1] + offset) corr = (*iter)[0] + 2;
            }

            Dihedral b(di()[0] + corr, di()[1] + offset, di()[2] + offset, di()[3] + offset);
            b.setType(di().type());

            // now we have to check whether we will add this one. We will not add
            // it if it 1) already exists or 2) refers to atoms less than nn
            int add = 1;

            // we do this by looking at the first three elements only, because
            // for a beginning group, we sometimes don't really know what the
            // fourth is.
            for (unsigned int j = 0; j < oldDihedral.size(); j++) {
                if (oldDihedral[j][0] == b[0] &&
                        oldDihedral[j][1] == b[1] &&
                        oldDihedral[j][2] == b[2]) {
                    // we found it, that means that it will not be added, but we have to
                    // change the fourth element in the one that we already had
                    oldDihedral[j][3] = b[3];
                    add = 0;
                }
            }
            // now check for the numbers in b
            // is this check really needed?
            // Yes! if you do things like ALA COO- NH3+ ALA
            //      it goes wrong otherwise
            if (b[0] < nn || b[1] < nn || b[2] < nn || b[3] < nn) add = 0;

            // here we handle the special case for CYS1 starting as the first residue
            if (b[1] == -6 + offset && b[2] == -5 + offset && b[3] == -4 + offset) add =1;
            if (b[0] == -5 + offset && b[1] == -6 + offset) add =1;
            if (b[3] == -6 + offset) add=1;

            // and if we still want it, we add it
            if (add) lt.addDihedral(b);
        }
        // now we have to re-add the old ones
        for (unsigned int j = 0; j < oldDihedral.size(); j++)
            lt.addDihedral(oldDihedral[j]);
    }
    // LJ exceptions
    LJExceptionIterator lji(bb);
    for (; lji; ++lji) {
        LJException lj(lji()[0] + offset, lji()[1] + offset);
        lj.setType(lji().type());
        lj.indicate() = lji().indicate();
        lj.cond() = lji().cond();

        //check if it exists already
        int found = 0;
        if (rep > 0) {
            for (std::set<gcore::LJException>::const_iterator iter = lt.ljexceptions().begin();
                    iter != lt.ljexceptions().end(); ++iter) {
                if ((*iter)[0] == lj[0] && (*iter)[1] == lj[1]) found = 1;
            }
        }
        if (!found)
            lt.addLJException(lj);
    }
}

int addBegin(gcore::LinearTopology &lt,
        BbSolute bb, int resnum) {
    int na = lt.atoms().size();

    if (na == 0) {
        //we just add all the atoms
        for (int i = 0; i < bb.numAtoms(); i++) {
            lt.addAtom(bb.atom(i));
            lt.setResNum(i, resnum);
        }
    } else {
        //we have to adapt the exclusions
        for (int i = 0; i < bb.numAtoms(); i++) {
            lt.addAtom(bb.atom(i));

            Exclusion e;
            for (int j = 0; j < bb.atom(i).exclusion().size(); j++)
                e.insert(bb.atom(i).exclusion().atom(j) + na);
            lt.atoms()[na + i].setExclusion(e);
            lt.setResNum(na + i, resnum);
        }
    }
    return bb.rep();
}

void addEnd(gcore::LinearTopology &lt,
        BbSolute bb, int resnum) {
    int strt = lt.atoms().size() + bb.rep();


    // first we check for any atoms in bb with a iac == -1
    vector<int> search;
    for(int i=0; i< bb.numAtoms(); i++){
       if ( bb.atom(i).iac() < -2 ) {
                    throw gromos::Exception("addCovEnd",
                        "IAC <-1 found in " + bb.resName() + " \n");
       }


      //if(bb.atom(i).iac()<-1){
      if(bb.atom(i).iac() == -2 ){
	search.push_back(i);
        std::cerr << "# WARNING\n" << "# For atom " << bb.atom(i).name() << " in MTBUILDBLEND group " << bb.resName() << " only the CHARGE\n" << "# is transferred to the last atom with this name in the chain." << std::endl;
      }
    }
    //we completely replace the last rep atoms,
    //but not the ones we need to search for
    for (unsigned int i = 0; i<(-bb.rep()-search.size()); i++)
        lt.atoms().pop_back();


    // now we search for the atoms based on the name
    for(unsigned int i=0; i<search.size(); i++){

      int j=lt.atoms().size()-1;
      while(lt.atoms()[j].name()!=bb.atom(search[i]).name())
        j--;

      // the atoms get the iac, charge and mass from the one in the endgroup
      // but it keeps its own chargegroup code and exclusions

      //lt.atoms()[j].setIac(-(bb.atom(search[i]).iac()+1)-1);
      lt.atoms()[j].setCharge(bb.atom(search[i]).charge()+lt.atoms()[j].charge());
      //lt.atoms()[j].setMass(bb.atom(search[i]).mass());
    }

    //and we add our new atoms
    for (int i = 0; i < bb.numAtoms(); i++) {
        // if we already did it as a search atom we skip it now
        if(bb.atom(i).iac()>=-1){
          Exclusion e;
          for (int j = 0; j < bb.atom(i).exclusion().size(); j++)
              e.insert(bb.atom(i).exclusion().atom(j) + strt);

          lt.addAtom(bb.atom(i));
          lt.atoms()[strt + i].setExclusion(e);
          lt.setResNum(strt + i, resnum);
      }
    }
}

void addCovEnd(gcore::LinearTopology &lt,
        BbSolute bb, int offset) {
    int found = 0;
    BondIterator bi(bb);

    for (; bi; ++bi) {
        Bond b(bi()[0] + offset, bi()[1] + offset);
        b.setType(bi().type());

        // if we are at the end of the tail
        if (bb.rep() < 0) {

            // first we see if the bond has negative values
            // In that case it should be present already
            // So in this case it might be a bit weird to search over the bonds for a bond.
            //if (bi()[0] < 0 || bb.atom(bi()[0]).iac()<-1) {
            if (bi()[0] < 0 || bb.atom(bi()[0]).iac() == -2 ) {
                std::set<int> candidates, atoms;
                atoms.insert(bi()[1] + offset);
                candidates = bondedAtoms(lt.bonds(), atoms, offset);

                if (candidates.size() == 0)
                    throw gromos::Exception("addCovEnd",
                        "Bond to the previous residue in " + bb.resName() + " is not found\n");
                if (candidates.size() != 1)
                    throw gromos::Exception("addCovEnd",
                        "Bond to the previous residue in " + bb.resName() + " is ambiguous\n");
                b[0] = *candidates.begin();
            }

            //search if this bond is present already
            found = 0;
            std::set<gcore::Bond>::iterator to_erase;
            for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
                    iter != lt.bonds().end(); ++iter) {
                if ((*iter)[0] == b[0] && (*iter)[1] == b[1]) {
                    to_erase = iter;
                    found = 1;
                }
            }

            if (found) lt.bonds().erase(to_erase);
            lt.addBond(b);
        }            // at the beginning of the tail, there is none of this crap
        else lt.addBond(b);
    }

    //now, the angles
    AngleIterator ai(bb);

    for (; ai; ++ai) {
        Angle b(ai()[0] + offset, ai()[1] + offset, ai()[2] + offset);
        b.setType(ai().type());

        //in case we are at the end of the chain, we should check for negative
        //values, it is still only the first that could be negative

        //if (bb.rep() < 0 && (ai()[0] < 0 || bb.atom(ai()[0]).iac() < -1)) {
        if (bb.rep() < 0 && (ai()[0] < 0 || bb.atom(ai()[0]).iac() ==  -2)) {
            std::set<int> candidates, atoms;
            atoms.insert(ai()[1] + offset);
            candidates = bondedAtoms(lt.bonds(), atoms, offset);

            if (candidates.size() == 0)
                throw gromos::Exception("addCovEnd",
                    "Angle to the previous residue in " + bb.resName() + " is not found\n");
            if (candidates.size() != 1)
                throw gromos::Exception("addCovEnd",
                    "Angle to the previous residue in " + bb.resName() + " is ambiguous\n");
            b[0] = *candidates.begin();
        }

        //search if this angle is present already
        std::vector<std::set<gcore::Angle>::iterator > to_erase;
        for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
                iter != lt.angles().end(); ++iter) {
            if ((*iter)[0] == b[0] &&
                    (*iter)[1] == b[1] &&
                    (*iter)[2] == b[2]) {
                to_erase.push_back(iter);
            }
        }
        for (unsigned int i = 0; i < to_erase.size(); ++i)
            lt.angles().erase(to_erase[i]);
        lt.addAngle(b);
    }

    //impropers
    ImproperIterator ii(bb);
    for (; ii; ++ii) {
        Improper b(ii()[0] + offset, ii()[1] + offset, ii()[2] + offset, ii()[3] + offset);
        b.setType(ii().type());

        // in case we are at the end of the chain, we should check for negative values
        // now it can be any one of the elements
        if (bb.rep() < 0) {
            std::set<int> candidates, atoms, negs;
            int list[4];

            for (int i = 0; i < 4; i++) {
                list[i] = ii()[i] + offset;
                if (ii()[i] >= 0 && bb.atom(ii()[i]).iac()>=-1 )
                    atoms.insert(ii()[i] + offset);
                else
                    negs.insert(i);
            }

            if (atoms.size() < 4) {
                candidates = bondedAtoms(lt.bonds(), atoms, offset);

                if (candidates.size() == 0)
                    throw gromos::Exception("addCovEnd",
                        "Improper to the previous residue in " + bb.resName() + " is not found\n");
                if (candidates.size() != 1)
                    throw gromos::Exception("addCovEnd",
                        "Improper to the previous residue in " + bb.resName() + " is ambiguous\n");
                if (negs.size() > 1)
                    throw gromos::Exception("addCovEnd",
                        "I'm not sure I can handle multiple negative values in " + bb.resName() + " impropers\n");
                // now set b to a new improper, because the ordering might have changed
                list[*negs.begin()] = *candidates.begin();
                b = Improper(list[0], list[1], list[2], list[3]);
                b.setType(ii().type());
            }
        }

        // search if this improper is present already, because of the 'random'
        // order in impropers, we will have to put them all in a set
        // and see if all elements are present
        std::vector<std::set<gcore::Improper>::iterator > to_erase;
        for (std::set<gcore::Improper>::iterator iter = lt.impropers().begin();
                iter != lt.impropers().end(); ++iter) {
            std::set<int> tryset;
            for (int j = 0; j < 4; j++) tryset.insert(b[j]);
            if (tryset.count((*iter)[0]) &&
                    tryset.count((*iter)[1]) &&
                    tryset.count((*iter)[2]) &&
                    tryset.count((*iter)[3])) {
                to_erase.push_back(iter);
            }
        }
        for (unsigned int i = 0; i < to_erase.size(); ++i)
            lt.impropers().erase(to_erase[i]);
        lt.addImproper(b);
    }

    //Dihedrals
    // the dihedrals should only be replaced from the original topology.
    // if the end group wants to add two dihedrals that are the same,
    // this should be allowed. So we have to keep track of the dihedrals
    // that are being added by this routine already.
    set<Dihedral> added_dihedrals;
    DihedralIterator di(bb);
    for (; di; ++di) {
        Dihedral b(di()[0] + offset, di()[1] + offset, di()[2] + offset, di()[3] + offset);
        b.setType(di().type());

        // in case we are at the end of the chain, we should check for negative values
        // now it can be any one of the elements
        if (bb.rep() < 0) {
            std::set<int> candidates, atoms, negs;
            int list[4];

            for (int i = 0; i < 4; i++) {
                list[i] = di()[i] + offset;
                if (di()[i] >= 0 && bb.atom(di()[i]).iac() >= -1) {
                    atoms.insert(di()[i] + offset);
                }
                else
                    negs.insert(i);
            }

            if (atoms.size() < 4) {

                candidates = bondedAtoms(lt.bonds(), atoms, offset);

                if (candidates.size() == 0)
                    throw gromos::Exception("addCovEnd",
                        "Dihedral to the previous residue in " + bb.resName() + " is not found\n");
                if (candidates.size() != 1)
                    std::cerr << "Warning: Dihedral to the previous residue in "
                        << bb.resName() << " is ambiguous" << std::endl;
                if (negs.size() > 1)
                    throw gromos::Exception("addCovEnd",
                        "I'm not sure I can handle multiple negative values in " + bb.resName() + " dihedrals\n");
                // now set b to a new dihedral, because the ordering might have changed
                list[*negs.begin()] = *candidates.begin();
                b = Dihedral(list[0], list[1], list[2], list[3]);
                b.setType(di().type());
            }
        }

        //search if this dihedral is present alread
        std::vector<std::set<gcore::Dihedral>::iterator > to_erase;
        for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
                iter != lt.dihedrals().end(); ++iter) {
            if ((*iter)[0] == b[0] &&
                    (*iter)[1] == b[1] &&
                    (*iter)[2] == b[2] &&
                    (*iter)[3] == b[3]) {
                // if it is not one that we added ourselves before
                if (added_dihedrals.count(*iter) == 0) {
                    to_erase.push_back(iter);
                }
            }
        }
        for (unsigned int i = 0; i < to_erase.size(); ++i)
            lt.dihedrals().erase(to_erase[i]);

        lt.addDihedral(b);
        added_dihedrals.insert(b);
    }
    // LJ Exceptions
    LJExceptionIterator lji(bb);
    for (; lji; ++lji) {
        LJException lj(lji()[0] + offset, lji()[1] + offset);
        lj.setType(lji().type());
        lj.indicate() = lji().indicate();
        lj.cond() = lji().cond();

        lt.addLJException(lj);
    }
}

void setCysteines(gcore::LinearTopology &lt,
        int a, int b) {
    // brace yourselves, this will be ugly!
    //
    // Make_top requires that
    // 1. both residues to be linked as Cysteine bridge have a CA
    //    (checked in make_top.cc)
    // 2. the atom that is being linked comes two atoms after CA
    // 3. the first linked atom of the second residue is referred to in the first
    //    residue with a number 8 lower than the atom number of the first CA
    //
    // exclusions
    for (int i = a + 1; i < a + 3; i++) {
        Exclusion e;
        for (int j = 0; j < lt.atoms()[i].exclusion().size(); j++) {
            if (lt.atoms()[i].exclusion().atom(j) < a)
                e.insert(a + b - 6 - lt.atoms()[i].exclusion().atom(j));
            else
                e.insert(lt.atoms()[i].exclusion().atom(j));
        }
        lt.atoms()[i].setExclusion(e);
    }

    // bond
    int added = 0, removed = 0, bt = 0;

    // we'll do this in two parts. First remove the bond to a negative number
    // and then add a new bond at an appropriate place (only if it was removed)
    for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
            iter != lt.bonds().end(); ++iter) {

        if ((*iter)[0] == a - 8 && (*iter)[1] == a + 2) {
            bt = iter->type();
            lt.bonds().erase(iter);
            removed = 1;
            break;
        }
    }
    if (!removed)
        throw gromos::Exception("setCysteines",
            "Bond connecting the building blocks was not found.\n\n"
            "Make_top requires that\n"
            "1. both residues to be linked as Cysteine bridge have a CA\n"
            "2. the atom that is being linked comes two atoms after CA\n"
            "3. the first linked atom of the second residue is referred to in\n"
            "   the first residue with a number 8 lower than the atom number of \n"
            "   the first CA");

    // wait with adding the corrected bond, to prevent later errors with
    // parsing
    for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
            iter != lt.bonds().end(); ++iter) {
        // we try to put it at the appropriate place
        if (!added && removed && (*iter)[0] == a + 1 && (*iter)[1] == a + 2) {
            added = 1;
            Bond bb(a + 2, b + 2);
            bb.setType(bt);
            lt.bonds().insert(iter, bb);
            break;
        }
    }

    //two angles
    int removed_a1 = 0, removed_a2 = 0;
    Angle bb(a + 1, a + 2, b + 2);
    for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
            iter != lt.angles().end(); ++iter) {
        if ((*iter)[0] == a - 8 && (*iter)[1] == a + 2 && (*iter)[2] == a + 1) {
            bb.setType(iter->type());
            lt.angles().erase(iter);
            removed_a1 = 1;
            break;
        }
    }
    if (removed_a1)
        lt.angles().insert(bb);
    else
        throw gromos::Exception("setCysteines", "Angle connecting the building blocks was not found");

    bb = Angle(a + 2, b + 2, b + 1);
    for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
            iter != lt.angles().end(); ++iter) {
        if ((*iter)[0] == a - 7 && (*iter)[1] == a - 8 && (*iter)[2] == a + 2) {
            bb.setType(iter->type());
            lt.angles().erase(iter);
            removed_a2 = 1;
            break;
        }
    }
    if (removed_a2)
        lt.angles().insert(bb);
    else
        throw gromos::Exception("setCysteines", "Angle connecting the building blocks was not found");

    //three dihedrals
    int removed_d1 = 0, removed_d2 = 0; // removed_d3 = 0;
    Dihedral di(a, a + 1, a + 2, b + 2);
    for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
            iter != lt.dihedrals().end(); ++iter) {
        if ((*iter)[0] == a && (*iter)[1] == a + 1 && (*iter)[2] == a + 2 && (*iter)[3] == a - 8) {
            di.setType(iter->type());
            lt.dihedrals().erase(iter);
            removed_d1 = 1;
            break;
        }
    }
    if (removed_d1)
        lt.dihedrals().insert(di);
    else
        throw gromos::Exception("setCysteines", "Connecting dihedral angle between two residues not found");

    di = Dihedral(a + 1, a + 2, b + 2, b + 1);
    for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
            iter != lt.dihedrals().end(); ++iter) {
        if ((*iter)[0] == a - 7 && (*iter)[1] == a - 8 && (*iter)[2] == a + 2 && (*iter)[3] == a + 1) {
            di.setType(iter->type());
            lt.dihedrals().erase(iter);
            removed_d2 = 1;
            break;
        }
    }
    if (removed_d2)
        lt.dihedrals().insert(di);
    else
        throw gromos::Exception("setCysteines", "Connecting dihedral angle between two residues not found");

    di = Dihedral(b, b + 1, b + 2, a + 2 );
    for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
            iter != lt.dihedrals().end(); ++iter) {
        if ((*iter)[0] == a + 2 && (*iter)[1] == a - 8 && (*iter)[2] == a - 7 && (*iter)[3] == a - 6) {
            di.setType(iter->type());
            lt.dihedrals().erase(iter);
            //removed_d3 = 1;
            break;
        }
    }
    if (removed_d2)
        lt.dihedrals().insert(di);
    else
        throw gromos::Exception("setCysteines", "Connecting dihedral angle between two residues not found");


}

void setHeme(gcore::LinearTopology &lt,
        int a1, int a2, int b) {
    // brace yourselves, this will be ugly!
    // we assume the following:
    // a1: atomnumber of the his1 - ca
    // a2: atomnumber of the his1 - ne2
    // b:  atomnumber of the first heme atom (fe)
    // exclusions
    for (int i = a1; i <= a2; i++) {
        Exclusion e;
        for (int j = 0; j < lt.atoms()[i].exclusion().size(); j++) {
            if (lt.atoms()[i].exclusion().atom(j) < a1)
                e.insert(b + a1 - 4 - lt.atoms()[i].exclusion().atom(j));
            else
                e.insert(lt.atoms()[i].exclusion().atom(j));
        }
        lt.atoms()[i].setExclusion(e);
    }
    // bonds
    int added = 0;
    vector<Bond> bonds_to_add;
    vector<std::set<gcore::Bond>::iterator> bonds_to_remove;
    for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
            iter != lt.bonds().end(); ++iter) {

        if ((*iter)[0] < a1 && (*iter)[1] == a2) {
            Bond bb(a2, b + a1 - 4 - (*iter)[0]);
            bb.setType(iter->type());
	    bonds_to_remove.push_back(iter);
            bonds_to_add.push_back(bb);
        }
    }
    for (unsigned int j = 0; j < bonds_to_remove.size(); j++)
      lt.bonds().erase(bonds_to_remove[j]);
    for (unsigned int j = 0; j < bonds_to_add.size(); j++)
      lt.bonds().insert(bonds_to_add[j]);
    //two kinds of angles
    vector<Angle> angles_to_add;
    vector<std::set<gcore::Angle>::iterator> angles_to_remove;
    for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
            iter != lt.angles().end(); ++iter) {

        if ((*iter)[0] < a1 && (*iter)[1] == a2) {
            Angle bb((*iter)[2], a2, b + a1 - 4 - (*iter)[0]);
            bb.setType(iter->type());
	    angles_to_remove.push_back(iter);
            angles_to_add.push_back(bb);
        }
        if ((*iter)[0] < a1 && (*iter)[1] < a1 && (*iter)[2] == a2) {
            Angle bb(a2, b + a1 - 4 - (*iter)[1], b + a1 - 4 - (*iter)[0]);
            bb.setType(iter->type());
	    angles_to_remove.push_back(iter);
            angles_to_add.push_back(bb);
        }
    }
    for (unsigned int j = 0; j < angles_to_remove.size(); j++)
        lt.angles().erase(angles_to_remove[j]);
    for (unsigned int j = 0; j < angles_to_add.size(); j++)
        lt.angles().insert(angles_to_add[j]);
    //no impropers?
    // dihedrals because it messes with the peptide linking, we had to
    // hardcode any found dihedrals to go to b, b+1
    // after peptide linking we have changed the b+1 number and it will
    // be virtually impossible to get back the original number
    vector<Dihedral> dihedrals_to_add;
    vector<std::set<gcore::Dihedral>::iterator> dihedrals_to_remove;
    for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
            iter != lt.dihedrals().end(); ++iter) {

        if ((*iter)[0] < a1 && (*iter)[1] < a1 && (*iter)[2] == a2) {
	    Dihedral bb((*iter)[3], a2, b, b+1);
            bb.setType(iter->type());
	    dihedrals_to_remove.push_back(iter);
	    dihedrals_to_add.push_back(bb);
        }
    }
    for (unsigned int j = 0; j < dihedrals_to_remove.size(); j++)
        lt.dihedrals().erase(dihedrals_to_remove[j]);
    for (unsigned int j = 0; j < dihedrals_to_add.size(); j++)
        lt.dihedrals().insert(dihedrals_to_add[j]);
}

std::set<int> bondedAtoms(std::set<gcore::Bond> &bonds,
        std::set<int> atoms,
        int offset) {
    std::set<int> candidates;
    int min = offset;
    //loop over the atoms
    for (std::set<int>::const_iterator iter = atoms.begin(), to = atoms.end();
            iter != to; ++iter)
        if (*iter < min)min = *iter;

    //loop over the atoms
    for (std::set<int>::const_iterator iter = atoms.begin(), to = atoms.end();
            iter != to; ++iter)
        //loop over the bonds
        for (std::set<gcore::Bond>::const_iterator ib = bonds.begin();
                ib != bonds.end(); ++ib)
            //if ((*ib)[1] == *iter && (*ib)[0] < min){
            if ((*ib)[1] == *iter && (*ib)[0] <= min) {
                // only add it to candidates if not in atoms set
                if ( atoms.count( ( (*ib)[0] ) )  == 0  ) {
                candidates.insert((*ib)[0]);
               }
            }
    return candidates;
}

void prepareCyclization(gcore::LinearTopology &lt) {
    // this function just adds three starting atoms to the lt,
    // connected by bonds.
    for (int i = 0; i < 3; i++) {
        AtomTopology at;
        at.setIac(18);
        lt.addAtom(at);
    }
    lt.addBond(Bond(0, 1));
    lt.addBond(Bond(1, 2));
    lt.addBond(Bond(1, 3));
}

void cyclize(gcore::LinearTopology &lt) {
    int na = lt.atoms().size();

    // first flag the first three atoms to be removed
    lt.atoms()[0].setIac(-1);
    lt.atoms()[1].setIac(-1);
    lt.atoms()[2].setIac(-1);
    // the last two atoms do not have exclusions yet
    Exclusion e_new;
    lt.atoms()[na - 2].setExclusion(e_new);
    lt.atoms()[na - 1].setExclusion(e_new);

    // the exclusions of atom 2 (na-1) have to be redistributed
    for (int i = 0; i < lt.atoms()[2].exclusion().size(); i++) {
        int excluded_atom = lt.atoms()[2].exclusion().atom(i);
        Exclusion e = lt.atoms()[excluded_atom].exclusion();
        e.insert(na - 1);
        lt.atoms()[excluded_atom].setExclusion(e);
    }

    // the exclusions of atom 1 (na-2) also have to be redistributed, but
    // if it is atom na-1 it is stored for atom na-2
    for (int i = 0; i < lt.atoms()[1].exclusion().size(); i++) {
        int excluded_atom = lt.atoms()[1].exclusion().atom(i);
        if (excluded_atom < 3) {
            Exclusion e;
            e.insert(excluded_atom + na - 3);
            lt.atoms()[na - 2].setExclusion(e);
        } else {
            Exclusion e = lt.atoms()[excluded_atom].exclusion();
            e.insert(na - 2);
            lt.atoms()[excluded_atom].setExclusion(e);
        }
    }

    // to do this nicely, we should loop over all remaining atoms and
    // redistribute all exclusions that are >= na
    // of course it is only relevant for the last few atoms.
    // we can make use of the fact that exclusions are sorted (are they still?)
    for (int i = 3; i < na; i++) {
        int n_excl = lt.atoms()[i].exclusion().size();
        if (n_excl && lt.atoms()[i].exclusion().atom(n_excl - 1) >= na) {
            Exclusion new_e_for_this_atom;
            for (int j = 0; j < n_excl; j++) {
                int excluded_atom = lt.atoms()[i].exclusion().atom(j);

                if (excluded_atom < na)
                    new_e_for_this_atom.insert(excluded_atom);
                else {
                    Exclusion new_e_for_ex_atom = lt.atoms()[excluded_atom - na + 3].exclusion();
                    new_e_for_ex_atom.insert(i);
                    lt.atoms()[excluded_atom - na + 3].setExclusion(new_e_for_ex_atom);
                }
            }
            lt.atoms()[i].setExclusion(new_e_for_this_atom);
        }
    }
    // find which atom maps to atom 0; this is the atom which is bonded to atom
    // na - 2 with a number less than na - 2
    int atomA = 0;

    for (std::set<gcore::Bond>::const_iterator iter = lt.bonds().begin();
            iter != lt.bonds().end(); ++iter) {
        if ((*iter)[1] == na - 2) atomA = (*iter)[0];
    }
    if (atomA == 0)
        throw (gromos::Exception("cyclise", "Cannot find atom A"));

    // covalent interactions!

    // loop over all bonds
    std::set<gcore::Bond> newBonds;
    for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
            iter != lt.bonds().end(); ++iter) {
        if ((*iter)[1] >= na) {
            //create a new bond (here don't complain about the order of atoms)
            Bond b((*iter)[0], (*iter)[1] - na + 3, 0);
            b.setType((*iter).type());
            newBonds.insert(b);
        } else {

            newBonds.insert(*iter);
        }

    }
    lt.bonds() = newBonds;

    // loop over all angles
    std::set<gcore::Angle> newAngles;
    for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
            iter != lt.angles().end(); ++iter) {
        if ((*iter)[2] >= na) {
            //create a new angle (here don't complain about the order of atoms)
            Angle a((*iter)[0], (*iter)[1], (*iter)[2] - na + 3, 0);
            a.setType((*iter).type());
            newAngles.insert(a);
        } else if ((*iter)[0] < 3) {
            //create a new angle (here don't complain about the order of atoms)
            Angle a((*iter)[0] + na - 3, (*iter)[1], (*iter)[2], 0);
            a.setType((*iter).type());
            newAngles.insert(a);
        } else newAngles.insert(*iter);
    }
    lt.angles() = newAngles;

    // loop over all impropers
    std::set<gcore::Improper> newImpropers;
    for (std::set<gcore::Improper>::iterator iter = lt.impropers().begin();
            iter != lt.impropers().end(); ++iter) {
        // anyone can be too high or too low
        int replace = 0;
        for (int i = 0; i < 4; i++) {
            if ((*iter)[i] < 3 || (*iter)[i] >= na) replace = 1;
        }
        int at[4];

        if (replace) {

            for (int i = 0; i < 4; i++) {
                if ((*iter)[i] < 3) at[i] = (*iter)[i] + na - 3;
                else if ((*iter)[i] >= na) at[i] = (*iter)[i] - na + 3;
                else at[i] = (*iter)[i];
            }
            Improper ii(at[0], at[1], at[2], at[3], 0);
            ii.setType(iter->type());
            newImpropers.insert(ii);
        } else {
            newImpropers.insert(*iter);
        }
    }
    lt.impropers() = newImpropers;

    // loop over all dihedrals
    std::set<gcore::Dihedral> newDihedrals;
    for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
            iter != lt.dihedrals().end(); ++iter) {
        // anyone can be too high or too low
        int replace = 0;
        for (int i = 0; i < 4; i++) {
            if ((*iter)[i] < 3 || (*iter)[i] >= na) replace = 1;
        }
        int at[4];

        if (replace) {

            for (int i = 0; i < 4; i++) {
                if ((*iter)[i] == 0) at[i] = atomA;
                else if ((*iter)[i] < 3) at[i] = (*iter)[i] + na - 3;
                else if ((*iter)[i] >= na) at[i] = (*iter)[i] - na + 3;
                else at[i] = (*iter)[i];
            }
            Dihedral d(at[0], at[1], at[2], at[3], 0);
            d.setType(iter->type());
            newDihedrals.insert(d);
        } else {
            newDihedrals.insert(*iter);
        }
    }
    lt.dihedrals() = newDihedrals;

    lt.removeAtoms();

    // Now check if the last charge group is closed.
    // If not, move atoms from top to bottom of list
    // untill the last charge group is closed
    // Note that moved atoms will be transfered from the first
    // residue to the last one (a warning will be written out)
    na = lt.atoms().size();
    int nrtransf = 0;
    while (lt.atoms()[na - 1].chargeGroup() != 1) {

        //add the first atom to the end
        lt.addAtom(lt.atoms()[0]);
        na = lt.atoms().size();
        lt.setResNum(na - 1, lt.resMap()[na - 2]);
        nrtransf++;

        //mark atom 1 to be removed
        lt.atoms()[0].setIac(-1);

        // now change the exclusions

        //last atom has no exclusion so make an empty one

        Exclusion e_new;
        lt.atoms()[na - 1].setExclusion(e_new);

        //exclusions of old atom 1 is then distributed
        // on the others (now as last exclusion)

        for (int i = 0; i < lt.atoms()[0].exclusion().size(); i++) {
            int excluded_atom = lt.atoms()[0].exclusion().atom(i);
            Exclusion e = lt.atoms()[excluded_atom].exclusion();
            e.insert(na - 1);
            lt.atoms()[excluded_atom].setExclusion(e);
        }

        //now loop through all bonds

        std::set<gcore::Bond> newBonds;
        for (std::set<gcore::Bond>::iterator iter = lt.bonds().begin();
                iter != lt.bonds().end(); ++iter) {
            if ((*iter)[0] == 0) {
                //create a new bond
                Bond b(na - 1, (*iter)[1]);
                b.setType((*iter).type());
                newBonds.insert(b);
            } else if ((*iter)[1] == 0) {
                //create a new bond
                Bond b((*iter)[0], na - 1);
                b.setType((*iter).type());
                newBonds.insert(b);
            } else {
                newBonds.insert(*iter);
            }
        }
        lt.bonds() = newBonds;

        //now loop through all angles

        std::set<gcore::Angle> newAngles;

        for (std::set<gcore::Angle>::iterator iter = lt.angles().begin();
                iter != lt.angles().end(); ++iter) {
            if ((*iter)[0] == 0) {
                //create a new angle
                Angle a(na - 1, (*iter)[1], (*iter)[2]);
                a.setType((*iter).type());
                newAngles.insert(a);
            }
            else if ((*iter)[1] == 0) {
                //create a new angle
                Angle a((*iter)[0], na - 1, (*iter)[2]);
                a.setType((*iter).type());
                newAngles.insert(a);
            } else if ((*iter)[2] == 0) {
                //create a new angle
                Angle a((*iter)[0], (*iter)[1], na - 1);
                a.setType((*iter).type());
                newAngles.insert(a);
            } else {
                newAngles.insert(*iter);
            }
        }
        lt.angles() = newAngles;


        //now loop through all improper dihedrals
        std::set<gcore::Improper> newImpropers;
        for (std::set<gcore::Improper>::iterator iter = lt.impropers().begin();
                iter != lt.impropers().end(); ++iter) {
            // anyone can be too high or too low
            int replace = 0;
            for (int i = 0; i < 4; i++) {
                if ((*iter)[i] == 0) replace = 1;
            }
            int at[4];

            if (replace) {
                for (int i = 0; i < 4; i++) {
                    if ((*iter)[i] == 0) {
                        at[i] = na - 1;
                    } else {
                        at[i] = (*iter)[i];
                    }
                }
                Improper ii(at[0], at[1], at[2], at[3]);
                ii.setType(iter->type());
                newImpropers.insert(ii);
            } else {
                newImpropers.insert(*iter);
            }
        }
        lt.impropers() = newImpropers;


        //now loop through all dihedrals
        std::set<gcore::Dihedral> newDihedrals;
        for (std::set<gcore::Dihedral>::iterator iter = lt.dihedrals().begin();
                iter != lt.dihedrals().end(); ++iter) {
            // anyone can be too high or too low
            int replace = 0;
            for (int i = 0; i < 4; i++) {
                if ((*iter)[i] == 0) replace = 1;
            }
            int at[4];

            if (replace) {

                for (int i = 0; i < 4; i++) {
                    if ((*iter)[i] == 0) {
                        at[i] = na - 1;
                    } else {
                        at[i] = (*iter)[i];
                    }
                }
                Dihedral d(at[0], at[1], at[2], at[3]);
                d.setType(iter->type());
                newDihedrals.insert(d);
            } else {
                newDihedrals.insert(*iter);
            }
        }
        lt.dihedrals() = newDihedrals;


        //finally remove the first atom

        lt.removeAtoms();
        na = lt.atoms().size();

    }
    if (nrtransf > 0) {
        std::cerr << "WARNING: in make_top::cyclize, ";
        std::cerr << nrtransf << " atom(s)" << std::endl;
        std::cerr << "were transfered from first residue to last one!" << std::endl;
    }

}
#endif
