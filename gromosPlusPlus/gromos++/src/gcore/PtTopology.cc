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

#include "PtTopology.h"


//gcore_PtTopology.cc
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include "AtomTopology.h"
#include "Molecule.h"
#include "LinearTopology.h"
#include "Bond.h"
#include "Angle.h"
#include "Improper.h"
#include "Dihedral.h"
#include "CrossDihedral.h"
#include "LJException.h"
#include "System.h"
#include "MoleculeTopology.h"
#include "Exclusion.h"
#include "../gromos/Exception.h"


namespace gcore
{
  
  PtTopology::PtTopology(const PtTopology & pt) {
    d_atomnum = pt.d_atomnum;
    d_residuenum = pt.d_residuenum;
    d_atomnames = pt.d_atomnames;
    d_residuenames = pt.d_residuenames;
    d_pertnames = pt.d_pertnames;
    d_iac = pt.d_iac;
    d_charge = pt.d_charge;
    d_mass = pt.d_mass;
    d_polarisability = pt.d_polarisability;
    d_dampingLevel = pt.d_dampingLevel;
    d_alphaLJ = pt.d_alphaLJ;
    d_alphaCRF = pt.d_alphaCRF;
    d_hasPolarisationParams = pt.d_hasPolarisationParams;
    d_multipt = pt.d_multipt;
    d_atompairs = pt.d_atompairs;
    d_bonds = pt.d_bonds;
    d_angles = pt.d_angles;
    d_impropers = pt.d_impropers;
    d_dihedrals = pt.d_dihedrals;
    d_crossdihedrals = pt.d_crossdihedrals;
  }
  
  void PtTopology::setSize(int a, int p)
  {
    d_iac.resize(p);
    d_mass.resize(p);
    d_charge.resize(p);
    d_pertnames.resize(p, "STATE");
    d_polarisability.resize(p);
    d_dampingLevel.resize(p);
    d_atompairs.resize(p);
    d_bonds.resize(p);
    d_angles.resize(p);
    d_impropers.resize(p);
    d_dihedrals.resize(p);
    d_crossdihedrals.resize(p);
    
    for(int i=0;i<p;i++){
      d_iac[i].resize(a, 0);
      d_mass[i].resize(a, 0.0);
      d_charge[i].resize(a, 0.0);
      d_polarisability[i].resize(a, 0.0);
      d_dampingLevel[i].resize(a, 0.0);
    }
    
    d_atomnames.resize(a, "AT");
    d_residuenames.resize(a, "RES");
    d_atomnum.resize(a, 0);
    d_residuenum.resize(a, 0);
    d_alphaLJ.resize(a, 0.0);
    d_alphaCRF.resize(a, 0.0);
  }
  
  void checkAtom(LinearTopology & top, unsigned int atom) {
    if (atom >= top.atoms().size()) {
      std::ostringstream os;
      os << "Atom out of range: " << atom + 1 << ".";
      throw gromos::Exception("PtTopology", os.str());
    }
  }
  
  void PtTopology::apply(System &sys, int iipt, int first) const
  {    
    LinearTopology top(sys);
    
    for (int i = 0; i < numAtoms(); ++i) {
      checkAtom(top, atomNum(i) + first);
      AtomTopology & atom = top.atoms()[atomNum(i) + first];
      if (atomName(i) != atom.name()) {
        std::cerr << "Warning: Atom names in (perturbation) topologies do not match" << std::endl
                << "Topology: " << atom.name()
                << "\tPerturbation topology: " << atomName(i)
                << "." << std::endl;
      }

      atom.setIac(iac(i, iipt));
      atom.setCharge(charge(i, iipt));
      
      if (d_multipt) {
        continue;
      }
      
      atom.setMass(mass(i, iipt));
      if (atom.isPolarisable() && hasPolarisationParameters()) {
        atom.setPolarisability(polarisability(i, iipt));
        atom.setDampingLevel(dampingLevel(i, iipt));
      }
      // loop to next atom in the perturbation list
    } // for atoms
    
    if (!d_multipt) {
      // do bonded and exclusions

      // apply the atom pairs
      for (std::set<AtomPairParam>::const_iterator it = atompairs(iipt).begin(),
              to = atompairs(iipt).end(); it != to; ++it) {
        const AtomPairParam & ap = *it;
        checkAtom(top, ap[0] + first);
        checkAtom(top, ap[1] + first);

        AtomTopology & atom_i = top.atoms()[ap[0] + first];
        switch (ap.param()) {
          case 0:
          { // excluded
            atom_i.exclusion().insert(ap[1] + first);
            atom_i.exclusion14().erase(ap[1] + first);
            break;
          }
          case 1:
          { // normal interaction
            atom_i.exclusion().erase(ap[1] + first);
            atom_i.exclusion14().erase(ap[1] + first);
            break;
          }
          case 2:
          { // 1-4 interaction
            atom_i.exclusion().erase(ap[1] + first);
            atom_i.exclusion14().insert(ap[1] + first);
            break;
          }
          default:
          {
            std::ostringstream os;
            os << "Invalid interaction type for atom pairs: 0..2 allowed but got "
                    << ap.param() << ".";
            throw gromos::Exception("PtTopology", os.str());
          }
        } // switch
      } // for atom pairs

      for (std::set<Bond>::const_iterator bond_it = bonds(iipt).begin(),
              bond_to = bonds(iipt).end(); bond_it != bond_to; ++bond_it) {
        Bond perturbed_bond((*bond_it)[0] + first, (*bond_it)[1] + first);
        perturbed_bond.setType(bond_it->type());
        checkAtom(top, perturbed_bond[0]);
        checkAtom(top, perturbed_bond[1]);

        // try to find the perturbed bond
        std::set<Bond>::iterator it = top.bonds().begin(), to = top.bonds().end();
        for (; it != to; ++it) {
          if ((*it)[0] == (perturbed_bond[0]) &&
              (*it)[1] == (perturbed_bond[1]))
            break;
        }
        if (it == to) {
          std::ostringstream os;
          os << "Bond " << perturbed_bond[0] + 1 << "-" 
             << perturbed_bond[1] + 1 
             << " does not exist.";
          throw gromos::Exception("PtTopology", os.str());
        }

        top.bonds().erase(it);
        top.bonds().insert(perturbed_bond);
      }

      for (std::set<Angle>::const_iterator angle_it = angles(iipt).begin(),
              angle_to = angles(iipt).end(); angle_it != angle_to; ++angle_it) {
        Angle perturbed_angle((*angle_it)[0] + first,(*angle_it)[1] + first,
			      (*angle_it)[2] + first);
        perturbed_angle.setType(angle_it->type());
        checkAtom(top, perturbed_angle[0]);
        checkAtom(top, perturbed_angle[1]);
        checkAtom(top, perturbed_angle[2]);

        // try to find the perturbed bond
        std::set<Angle>::iterator it = top.angles().begin(), to = top.angles().end();
        for (; it != to; ++it) {
          if ((*it)[0] == (perturbed_angle[0]) &&
              (*it)[1] == (perturbed_angle[1]) &&
              (*it)[2] == (perturbed_angle[2]))
            break;
        }
        if (it == to) {
          std::ostringstream os;
          os << "Angle " << perturbed_angle[0] + 1 
	      << "-" << perturbed_angle[1] + 1 
              << "-" << perturbed_angle[2] + 1  << " does not exist.";
          throw gromos::Exception("PtTopology", os.str());
        }

        top.angles().erase(it);
        top.angles().insert(perturbed_angle);
      }

      for (std::set<Improper>::const_iterator improper_it = impropers(iipt).begin(),
              improper_to = impropers(iipt).end(); improper_it != improper_to; ++improper_it) {
        Improper perturbed_improper((*improper_it)[0] + first, 
                                    (*improper_it)[1] + first,
                                    (*improper_it)[2] + first,
                                    (*improper_it)[3] + first);
        perturbed_improper.setType(improper_it->type());
        checkAtom(top, perturbed_improper[0]);
        checkAtom(top, perturbed_improper[1]);
        checkAtom(top, perturbed_improper[2]);
        checkAtom(top, perturbed_improper[3]);

        // try to find the perturbed bond
        std::set<Improper>::iterator it = top.impropers().begin(), to = top.impropers().end();
        for (; it != to; ++it) {
          if ((*it)[0] == (perturbed_improper[0]) &&
              (*it)[1] == (perturbed_improper[1]) &&
              (*it)[2] == (perturbed_improper[2]) &&
              (*it)[3] == (perturbed_improper[3]))
            break;
        }
        if (it == to) {
          std::ostringstream os;
          os << "Improper " << perturbed_improper[0] + 1 << "-" 
             << perturbed_improper[1] + 1 << "-" << perturbed_improper[2] + 1
             << "-" << perturbed_improper[3] + 1 << " does not exist.";
          throw gromos::Exception("PtTopology", os.str());
        }

        top.impropers().erase(it);
        top.impropers().insert(perturbed_improper);
      }
      for (std::set<Dihedral>::const_iterator dihedral_it = dihedrals(iipt).begin(),
              dihedral_to = dihedrals(iipt).end(); dihedral_it != dihedral_to; ++dihedral_it) {
        Dihedral perturbed_dihedral((*dihedral_it)[0] + first,
                                    (*dihedral_it)[1] + first,
                                    (*dihedral_it)[2] + first,
                                    (*dihedral_it)[3] + first);
        perturbed_dihedral.setType(dihedral_it->type());
        checkAtom(top, perturbed_dihedral[0]);
        checkAtom(top, perturbed_dihedral[1]);
        checkAtom(top, perturbed_dihedral[2]);
        checkAtom(top, perturbed_dihedral[3]);

        // try to find the perturbed bond
        std::set<Dihedral>::iterator it = top.dihedrals().begin(), to = top.dihedrals().end();
        for (; it != to; ++it) {
          if ((*it)[0] == (perturbed_dihedral[0]) &&
              (*it)[1] == (perturbed_dihedral[1]) &&
              (*it)[2] == (perturbed_dihedral[2]) &&
              (*it)[3] == (perturbed_dihedral[3]))
            break;
        }
        if (it == to) {
          std::ostringstream os;
          os << "Dihedral " << perturbed_dihedral[0] + 1 << "-"
             << perturbed_dihedral[1] + 1 << "-" << perturbed_dihedral[2] + 1
             << "-" << perturbed_dihedral[3] + 1 << " does not exist.";
          throw gromos::Exception("PtTopology", os.str());
        }

        top.dihedrals().erase(it);
        top.dihedrals().insert(perturbed_dihedral);
      }
      for (std::set<CrossDihedral>::const_iterator crossdihedral_it =  crossdihedrals(iipt).begin(),
              crossdihedral_to = crossdihedrals(iipt).end(); crossdihedral_it != crossdihedral_to; ++crossdihedral_it) {
        CrossDihedral perturbed_crossdihedral((*crossdihedral_it)[0] + first,
                                              (*crossdihedral_it)[1] + first,
                                              (*crossdihedral_it)[2] + first,
                                              (*crossdihedral_it)[3] + first,
                                              (*crossdihedral_it)[4] + first,
                                              (*crossdihedral_it)[5] + first,
                                              (*crossdihedral_it)[6] + first,
                                              (*crossdihedral_it)[7] + first);
        perturbed_crossdihedral.setType(crossdihedral_it->type());
        for(unsigned int i = 0; i < 8; ++i)
          checkAtom(top, perturbed_crossdihedral[i] + first);

        // try to find the perturbed cross dihedral
        std::set<CrossDihedral>::iterator it = top.crossdihedrals().begin(), to = top.crossdihedrals().end();
        for (; it != to; ++it) {
          if ((*it)[0] == (perturbed_crossdihedral[0]) &&
              (*it)[1] == (perturbed_crossdihedral[1]) &&
              (*it)[2] == (perturbed_crossdihedral[2]) &&
              (*it)[3] == (perturbed_crossdihedral[3]) &&
              (*it)[4] == (perturbed_crossdihedral[4]) &&
              (*it)[5] == (perturbed_crossdihedral[5]) &&
              (*it)[6] == (perturbed_crossdihedral[6]) &&
              (*it)[7] == (perturbed_crossdihedral[7]))
            break;
        }
        if (it == to) {
          std::ostringstream os;
          os << "CrossDihedral " << perturbed_crossdihedral[0] + 1 << "-" 
                                 << perturbed_crossdihedral[1] + 1 << "-" 
                                 << perturbed_crossdihedral[2] + 1 << "-"
                                 << perturbed_crossdihedral[3] + 1 << "-" 
                                 << perturbed_crossdihedral[4] + 1 << "-"
                                 << perturbed_crossdihedral[5] + 1 << "-" 
                                 << perturbed_crossdihedral[6] + 1 << "-" 
                                 << perturbed_crossdihedral[7] + 1 
             << " does not exist.";
          throw gromos::Exception("PtTopology", os.str());
        }

        top.crossdihedrals().erase(it);
        top.crossdihedrals().insert(perturbed_crossdihedral);
      }
    } // do bonded
    
    // the topology is adopted now. Apply the changes to the system
    System new_sys(top.parse());
    // copy over the data
    for(int m = 0; m < sys.numMolecules(); ++m) {
      sys.mol(m).topology() = new_sys.mol(m).topology();
    }
  }
  
  void PtTopology::setIac(int a, int p, int iac)
  {
    d_iac[p][a]=iac;
  }

  void PtTopology::setCharge(int a, int p, double q)
  {
    d_charge[p][a]=q;
  }
  
  void PtTopology::setMass(int a, int p, double q)
  {
    d_mass[p][a]=q;
  }
  
  void PtTopology::setPolarisability(int a, int p, double q)
  {
    d_polarisability[p][a]=q;
  }
  
  void PtTopology::setDampingLevel(int a, int p, double q)
  {
    d_dampingLevel[p][a]=q;
  }
  
  void PtTopology::setAtomName(int a, std::string name)
  {
    d_atomnames[a]=name;
  }
  void PtTopology::setResidueName(int a, std::string name)
  {
    d_residuenames[a]=name;
  }
  void PtTopology::setPertName(int p, std::string name)
  {
    d_pertnames[p]=name;
  }

  void PtTopology::setAtomNum(int a, int num)
  {
    d_atomnum[a]=num;
  }
  
  void PtTopology::setResidueNum(int a, int num)
  {
    d_residuenum[a]=num;
  }
  
  void PtTopology::setAlphaLJ(int a, double alpha)
  {
    d_alphaLJ[a]=alpha;
  }
  
  void PtTopology::setAlphaCRF(int a, double alpha)
  {
    d_alphaCRF[a]=alpha;
  }
  
  void PtTopology::setHasPolarisationParameters(bool b)
  {
    d_hasPolarisationParams=b;
  }
  
  void PtTopology::setMultiPt(bool multipt)
  {
    d_multipt=multipt;
  }
  
  PtTopology::PtTopology(System & sysA, System & sysB, bool quiet) {
    d_hasPolarisationParams = false; // we can set it to true later
    d_multipt = false;
    LinearTopology topA(sysA), topB(sysB);
    if (topA.atoms().size() != topB.atoms().size())
      throw gromos::Exception("PtTopology", "Both topologies need to have "
              "the same number of atoms.");
    
    setSize(0, 2);
    
     // loop over all atoms
    for(unsigned int i = 0; i < topA.atoms().size(); i++) {
      AtomTopology& atomA = topA.atoms()[i], atomB = topB.atoms()[i];
      
      bool polarisabilityChange = false;
      if (atomA.isPolarisable() && !atomB.isPolarisable()) {
        if (!quiet)
          std::cerr << "Warning: atom " << i+1 << " A is polarizable but B is not." << std::endl;
        atomB.setPolarisable(true);
        atomB.setPolarisability(0.0);
        atomB.setDampingLevel(0.0);
        atomB.setCosCharge(atomA.cosCharge());
        atomB.setDampingPower(atomA.dampingPower());
        polarisabilityChange = true;
      } else if (!atomA.isPolarisable() && atomB.isPolarisable()) {
        if (!quiet)
          std::cerr << "Warning: atom " << i+1 << " B is polarizable but A is not." << std::endl;
        atomA.setPolarisable(true);
        atomA.setPolarisability(0.0);
        atomA.setDampingLevel(0.0);
        atomA.setCosCharge(atomB.cosCharge());
        atomA.setDampingPower(atomB.dampingPower());
        polarisabilityChange = true;
      } else if (atomA.isPolarisable() && atomB.isPolarisable()) {
        // both are polarisable and their parameters may have changed
        if (atomA.polarisability() != atomB.polarisability() ||
            atomA.dampingLevel() != atomB.dampingLevel())
          polarisabilityChange = true;
        else
          polarisabilityChange = false;
      }
      
      if (atomA.iac() != atomB.iac() ||
          atomA.mass() != atomB.mass() ||
          atomA.charge() != atomB.charge() || 
          polarisabilityChange) {
        if (atomA.name() != atomB.name()) {
          // this is not bad but worth a warning
          if (!quiet)
            std::cerr << "Warning: atom " << i+1 << " name mismatch (" <<
                  atomA.name() << " != " << atomB.name() << ")" << std::endl;        
        }
        
        const int index = numAtoms();
        setSize(index+1, 2);
        setAtomNum(index, i);
        setResidueNum(index, topA.resMap()[i]);
        setAtomName(index, atomA.name());
        setResidueName(index, topA.resNames()[topA.resMap()[i]]);
        setIac(index, 0, atomA.iac());
        setIac(index, 1, atomB.iac());
        setMass(index, 0, atomA.mass());
        setMass(index, 1, atomB.mass());
        setCharge(index, 0, atomA.charge());
        setCharge(index, 1, atomB.charge());
        d_hasPolarisationParams = (d_hasPolarisationParams || polarisabilityChange);
        setPolarisability(index, 0, atomA.polarisability());
        setPolarisability(index, 1, atomB.polarisability());
        setDampingLevel(index, 0, atomA.dampingLevel());
        setDampingLevel(index, 1, atomB.dampingLevel());
      } // if perturbed atom
      
      // exclusions are difficult:
      // check wheter all atoms excluded by a certain atom in topology A
      // are also excluded in topology B and vice versa. 
      // don't forget 1-4 exclusions.
      const Exclusion& exA = atomA.exclusion(), exB = atomB.exclusion(),
                 ex14A = atomA.exclusion14(), ex14B = atomB.exclusion14();
      // exclusions of atomA also in atomB?
      for(int j = 0; j < exA.size(); j++) {
        unsigned stateA = 0; // excluded
        unsigned stateB = 1; // normal interaction
        if (exB.contains(exA.atom(j)))
          stateB = 0;
        else if (ex14B.contains(exA.atom(j)))
          stateB = 2;
        
        if (stateA != stateB) {
          atompairs(0).insert(AtomPairParam(i, exA.atom(j), stateA));
          atompairs(1).insert(AtomPairParam(i, exA.atom(j), stateB));
        }                 
      }
      // exclusions14 of atomA also in atomB?
      for(int j = 0; j < ex14A.size(); j++) {
        unsigned stateA = 2; // 1-4 excluded
        unsigned stateB = 1; // normal interaction
        if (exB.contains(ex14A.atom(j)))
          stateB = 0;
        else if (ex14B.contains(ex14A.atom(j)))
          stateB = 2;
        
        if (stateA != stateB) {
          atompairs(0).insert(AtomPairParam(i, ex14A.atom(j), stateA));
          atompairs(1).insert(AtomPairParam(i, ex14A.atom(j), stateB));
        }                 
      }
      // and the same the other way round:
      // exclusions of atomB also in atomA?
     for(int j = 0; j < exB.size(); j++) {
        unsigned stateB = 0; // excluded
        unsigned stateA = 1; // normal interaction
        if (exA.contains(exB.atom(j)))
          stateA = 0;
        else if (ex14A.contains(exB.atom(j)))
          stateA = 2;
        
        if (stateA != stateB) {
          atompairs(0).insert(AtomPairParam(i, exB.atom(j), stateA));
          atompairs(1).insert(AtomPairParam(i, exB.atom(j), stateB));
        }                 
      }
      // exclusions14 of atomB also in atomA?
      // Pascal: This was a bug - size of ex14A instead of ex14B
      // for(int j = 0; j < ex14A.size(); j++) {
      for(int j = 0; j < ex14B.size(); j++) {
        unsigned stateB = 2; // 1-4 excluded
        unsigned stateA = 1; // normal interaction
        if (exA.contains(ex14B.atom(j)))
          stateA = 0;
        else if (ex14A.contains(ex14B.atom(j)))
          stateA = 2;
        
        if (stateA != stateB) {
          atompairs(0).insert(AtomPairParam(i, ex14B.atom(j), stateA));
          atompairs(1).insert(AtomPairParam(i, ex14B.atom(j), stateB));
        }                 
      }
    } // for atoms
    
    if (!quiet && topA.bonds().size() != topB.bonds().size()) {
      std::cerr << "Warning: the topologies don't have the same "
              "number of bonds!" << std::endl
           << "         Taking bonds of state A and searching "
              "for changes in state B." << std::endl;
    }
    
    // loop over all bonds
    std::set<Bond>::const_iterator itBondA = topA.bonds().begin(), 
            toBondA = topA.bonds().end();
    
    for(; itBondA != toBondA; ++itBondA) {
      const Bond& bondA = *itBondA;
      bool found = false;
      
      // search for bondA in topology B -> loop over bonds in topology B
      std::set<Bond>::const_iterator itBondB = topB.bonds().begin(), 
            toBondB = topB.bonds().end();
      
      for(; itBondB != toBondB; ++itBondB) {
        const Bond& bondB = *itBondB;
        if (bondA[0] == bondB[0] && bondA[1] == bondB[1]) {
          found = true;
          if (bondA.type() != bondB.type()) {
            // add it
            bonds(0).insert(bondA);
            bonds(1).insert(bondB);
          }
          break;
        }
      }
      if (!quiet && !found) {
        std::cerr << "Warning: could not find bondA " << bondA[0]+1 << "-"
             << bondA[1]+1 << std::endl;
      }
    } // for bonds
    
    if (!quiet && topA.angles().size() != topB.angles().size()) {
      std::cerr << "Warning: the topologies don't have the same "
              "number of angles!" << std::endl
           << "         Taking bond angles of state A and searching "
              "for changes in state B." << std::endl;
    }
    
    // loop over all angles
    std::set<Angle>::const_iterator itAngleA = topA.angles().begin(), 
            toAngleA = topA.angles().end();
    
    for(; itAngleA != toAngleA; ++itAngleA) {
      const Angle& angleA = *itAngleA;
      bool found = false;
      
      // search for angleA in topology B
      std::set<Angle>::const_iterator itAngleB = topB.angles().begin(), 
            toAngleB = topB.angles().end();
      
      for(; itAngleB != toAngleB; ++itAngleB) {
        const Angle& angleB = *itAngleB;
        if (angleA[0] == angleB[0] &&
            angleA[1] == angleB[1] &&
            angleA[2] == angleB[2]) {
          found = true;
          if (angleA.type() != angleB.type()) {
            angles(0).insert(angleA);
            angles(1).insert(angleB);
          }
          break;
        }
      }
      if (!quiet && !found) {
        std::cerr << "Warning: could not find angleA " << angleA[0]+1 << "-"
             << angleA[1]+1 << "-" << angleA[2]+1 << std::endl;
      }
    } // for angles
    
    if (!quiet && topA.impropers().size() != topB.impropers().size()) {
      std::cerr << "Warning: the topologies don't have the same "
              "number of improper dihedrals!" << std::endl
           << "         Taking improper dihedrals of state A and searching "
              "for changes in state B." << std::endl;
    }
    
    // loop over impropers
    std::set<Improper>::const_iterator itImproperA = topA.impropers().begin(), 
            toImproperA = topA.impropers().end();
    
    for(; itImproperA != toImproperA; ++itImproperA) {
      const Improper& improperA = *itImproperA;
      bool found = false;
      
      // search for improperA in topology B
      std::set<Improper>::const_iterator itImproperB = topB.impropers().begin(), 
            toImproperB = topB.impropers().end();
      
      for(; itImproperB != toImproperB; ++itImproperB) {
        const Improper& improperB = *itImproperB;
        if (improperA[0] == improperB[0] &&
            improperA[1] == improperB[1] &&
            improperA[2] == improperB[2] &&
            improperA[3] == improperB[3]) {
          found = true;
          if (improperA.type() != improperB.type()) {
            impropers(0).insert(improperA);
            impropers(1).insert(improperB);
          }
          break;
        }
      }
      if (!quiet && !found) {
        std::cerr << "Warning: could not find improperA " << improperA[0]+1 << "-"
             << improperA[1]+1 << "-" << improperA[2]+1 << std::endl;
      }
    } // for impropers
    
    if (topA.dihedrals().size() != topB.dihedrals().size()) {
      std::cerr << "Warning: the topologies don't have the same "
              "number of dihedral dihedrals!" << std::endl
           << "         Taking dihedrals of state A and searching "
              "for changes in state B." << std::endl;
    }
    
    // loop over dihedrals
    std::set<Dihedral>::const_iterator itDihedralA = topA.dihedrals().begin(), 
            toDihedralA = topA.dihedrals().end();
    
    for(; itDihedralA != toDihedralA; ++itDihedralA) {
      const Dihedral& dihedralA = *itDihedralA;
      bool found = false;
      
      // search for dihedralA in topology B
      
      if (topB.dihedrals().find(dihedralA) != topB.dihedrals().end())
        continue;
      
      std::set<Dihedral>::const_iterator itDihedralB = topB.dihedrals().begin(), 
            toDihedralB = topB.dihedrals().end();
      
      for(; itDihedralB != toDihedralB; ++itDihedralB) {
        const Dihedral& dihedralB = *itDihedralB;
        if (dihedralA[0] == dihedralB[0] &&
            dihedralA[1] == dihedralB[1] &&
            dihedralA[2] == dihedralB[2] &&
            dihedralA[3] == dihedralB[3]) {
          found = true;
          if (dihedralA.type() != dihedralB.type()) {
            dihedrals(0).insert(dihedralA);
            dihedrals(1).insert(dihedralB);
          }
          break;
        }
      }
      if (!quiet && !found) {
        std::cerr << "Warning: could not find dihedral " << dihedralA[0]+1 << "-"
             << dihedralA[1]+1 << "-" << dihedralA[2]+1 << "-" << dihedralA[3]+1 << std::endl;
      }
    } // for dihedrals

    if (topA.crossdihedrals().size() != topB.crossdihedrals().size()) {
      std::cerr << "Warning: the topologies don't have the same "
              "number of cross-dihedrals!" << std::endl
           << "         Taking cross dihedrals of state A and searching "
              "for changes in state B." << std::endl;
    }

    // loop over cross dihedrals
    {
      std::set<CrossDihedral>::const_iterator itDihedralA = topA.crossdihedrals().begin(),
              toDihedralA = topA.crossdihedrals().end();

      for (; itDihedralA != toDihedralA; ++itDihedralA) {
        const CrossDihedral& dihedralA = *itDihedralA;
        bool found = false;

        // search for dihedralA in topology B

        std::set<CrossDihedral>::const_iterator itDihedralB = topB.crossdihedrals().begin(),
                toDihedralB = topB.crossdihedrals().end();

        for (; itDihedralB != toDihedralB; ++itDihedralB) {
          const CrossDihedral& dihedralB = *itDihedralB;
          if (dihedralA[0] == dihedralB[0] &&
                  dihedralA[1] == dihedralB[1] &&
	      dihedralA[2] == dihedralB[2] &&
                  dihedralA[3] == dihedralB[3] &&
                  dihedralA[4] == dihedralB[4] &&
                  dihedralA[5] == dihedralB[5] &&
                  dihedralA[6] == dihedralB[6] &&
                  dihedralA[7] == dihedralB[7]) {
            found = true;
            if (dihedralA.type() != dihedralB.type()) {
              crossdihedrals(0).insert(dihedralA);
              crossdihedrals(1).insert(dihedralB);
            }
            break;
          }
        }
        if (!quiet && !found) {
          std::cerr << "Warning: could not find cross dihedral "
                  << dihedralA[0] + 1 << "-" << dihedralA[1] + 1 << "-" << dihedralA[2] + 1 << "-" << dihedralA[3] + 1
                  << dihedralA[4] + 1 << "-" << dihedralA[5] + 1 << "-" << dihedralA[6] + 1 << "-" << dihedralA[7] + 1 << std::endl;
        }
      }
    }
    
    // for LJ exceptions
    if (topA.ljexceptions().size() != topB.ljexceptions().size()) {
      throw gromos::Exception("PtTopology", "The topologies don't have the "
			      "same number of LJ exceptions!\n LJ Exceptions "
			      "cannot be perturbed!");
    }

    // loop over LJ exceptions
    {
      std::set<LJException>::const_iterator itLJExceptionA = topA.ljexceptions().begin(),
              toLJExceptionA = topA.ljexceptions().end();

      for (; itLJExceptionA != toLJExceptionA; ++itLJExceptionA) {
        const LJException& ljexceptionA = *itLJExceptionA;
        bool found = false;

        // search for ljexceptionA in topology B

        std::set<LJException>::const_iterator itLJExceptionB = topB.ljexceptions().begin(),
                toLJExceptionB = topB.ljexceptions().end();

        for (; itLJExceptionB != toLJExceptionB; ++itLJExceptionB) {
          const LJException& ljexceptionB = *itLJExceptionB;
          if (ljexceptionA[0] == ljexceptionB[0] &&
	      ljexceptionA[1] == ljexceptionB[1]) {
            found = true;
            if (ljexceptionA.type() != ljexceptionB.type()) {
              throw gromos::Exception("PtTopology", "Different types found for "
				      "LJ Exception!\n LJ Exceptions cannot be "
				      " perturbed");
            }
            break;
          }
        }
        if (!quiet && !found) {
          std::cerr << "Warning: could not find LJ Exception "
                  << ljexceptionA[0] + 1 << "-" << ljexceptionA[1] + 1 << std::endl;
        }
      } // for LJExceptions
    }
  }
  
  PtTopology::PtTopology(std::vector<System*> & sys, bool quiet) {
    if (sys.empty())
      throw gromos::Exception("PtTopology", "Give topologies!");
      
    d_hasPolarisationParams = false; 
    d_multipt = true;
    std::vector<LinearTopology*> top;
    
    for(unsigned int i = 0; i < sys.size(); ++i)
      top.push_back(new LinearTopology(*sys[i]));
    
    const unsigned int num_atoms = top[0]->atoms().size();
    for(unsigned int i = 1; i < sys.size(); ++i) {
      if (num_atoms != top[i]->atoms().size())
        throw gromos::Exception("PtTopology", "Topologies need to have "
                "the same number of atoms.");
    }
    
     // loop over all atoms and save the check
    std::vector<int> changed;
    for(unsigned int i = 0; i < num_atoms; i++) {
      AtomTopology& atomA = top[0]->atoms()[i];
      for(unsigned int j = 1; j < sys.size(); ++j) {
        AtomTopology& atomB = top[j]->atoms()[i];
        // only check charge and iac for multiple perturbation
        if (atomA.iac() != atomB.iac() ||
            atomA.charge() != atomB.charge()) {
          if (atomA.name() != atomB.name()) {
            // this is not bad but worth a warning
            if (!quiet)
               std::cerr << "Warning: atom " << i+1 << " name mismatch (" <<
                    atomA.name() << " != " << atomB.name() << ")" << std::endl;
          }
          // it has changed
          changed.push_back(i);
        } // mismatch
      } // for systems
    } // for atoms
    
    setSize(changed.size(), sys.size());
    for(unsigned int j = 0; j < sys.size(); ++j) {
      std::ostringstream name;
      name << "S" << std::setfill('0') << std::setw(3) << j+1;
      setPertName(j, name.str());
    }
    
    for(unsigned int index = 0; index < changed.size(); index++) {
      const int i = changed[index];
      // set the static information
      setAtomNum(index, i);
      setResidueNum(index, top[0]->resMap()[i]);
      setAtomName(index, top[0]->atoms()[i].name());
      setResidueName(index, top[0]->resNames()[top[0]->resMap()[i]]);
      // set the parameter for all states
      for(unsigned int j = 0; j < sys.size(); ++j) {
        AtomTopology& atom = top[j]->atoms()[i];
        setIac(index, j, atom.iac());
        setCharge(index, j, atom.charge());
      } // for systems
    } // for atoms
  }

}



