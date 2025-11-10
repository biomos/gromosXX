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

// gcore_LinearTopology.h

#ifndef INCLUDED_GCORE_LINEARTOPOLOGY
#define INCLUDED_GCORE_LINEARTOPOLOGY

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Constraint.h"
#include "LJException.h"
#include "Improper.h"
#include "Bond.h"
#include "Angle.h"
#include "AtomTopology.h"
#include "CrossDihedral.h"
#include "System.h"
#include "Dihedral.h"

namespace gcore {

  class System;
  class AtomTopology;
  class Bond;
  class Constraint;
  class Angle;
  class Dihedral;
  class CrossDihedral;
  class Improper;
  //class LJException;

  /**
   * Class LinearTopology
   * Purpose: A preliminary container for the topology of a system
   *
   * Description:
   * A linearized topology, which does not know about molecules or solvent
   * but is useful when building up a topology, or changing a topology such 
   * that the Molecules change.
   *
   * @class LinearTopology
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::AtomTopology
   * @sa gcore::Bond
   * @sa gcore::Constraint
   * @sa gcore::Angle
   * @sa gcore::Improper
   * @sa gcore::Dihedral
   * @sa gcore::CrossDihedral
   */
  class LinearTopology {
    std::vector<AtomTopology> d_atom;
    std::set<Bond> d_bond;
    std::set<Constraint> d_constraints;
    std::set<Bond> d_dipole_bond;
    std::set<Angle> d_angle;
    std::set<Dihedral> d_dihedral;
    std::set<CrossDihedral> d_crossdihedral;
    std::set<Improper> d_improper;
    std::set<LJException> d_ljexception;
    std::vector<std::string> d_resname;
    std::map<int, int> d_resmap;

  public:
    /**
     * LinearTopology constructor, creates an empty LinearTopology
     */
    LinearTopology();
    /**
     * LinearTopology constructor, linearizes the topological information
     * of the system
     */
    LinearTopology(gcore::System &sys);
    /**
     * LinearTopology deconstructor
     */
    ~LinearTopology();

    // Methods
    /**
     * Method to add an atom to the MoleculeTopology
     * @param a An AtomTopology that is to be added; should be complete
     *          already
     */
    void addAtom(const AtomTopology &a);
    /**
     * Method to add a Bond to the MoleculeTopology
     * @param b The Bond that is to be added; should be complete already
     */
    void addBond(const Bond &b);
    /**
     * Method to add a Constraint to the MoleculeTopology
     * @param b The Constraint that is to be added; should be complete already
     */
    void addConstraint(const Constraint &b);
    /**
     * Method to add a Dipole Bond to the MoleculeTopology
     * @param b The Bond that is to be added; should be complete already
     */
    void addDipoleBond(const Bond &b);    
    /**
     * Method to add an Angle to the MoleculeTopology
     * @param b The Angle that is to be added; should be complete already
     */
    void addAngle(const Angle &b);
    /**
     * Method to add a Dihedral to the MoleculeTopology
     * @param b The Dihedral that is to be added; should be complete already
     */
    void addDihedral(const Dihedral &b);
    /**
     * Method to add a CrossDihedral to the MoleculeTopology
     * @param b The CrossDihedral that is to be added; should be complete already
     */
    void addCrossDihedral(const CrossDihedral &b);
    /**
     * Method to add an Improper to the MoleculeTopology
     * @param b The Improper that is to be added; should be complete already
     */
    void addImproper(const Improper &b);
    /**
     * Method to add a LJ exception to the MoleculeTopology
     * @param lj The LJ exception that is to be added; should be complete already
     */
    void addLJException(const LJException &lj);
    /**
     * Method to set the residue name
     * 
     * Within the MoleculeTopology a list of residues is kept so that you
     * can assign specific atoms to specific residues
     * @param res The number of the residue that gets a name
     * @param s   The name of this residue
     */
    void setResName(int res, const std::string &s);
    /**
     * Method to assign an atom to a residue
     *
     * Within the MoleculeTopology a list of residues is known, this method
     * allows you to assign an atom to a specific residue
     * @param atom The atom number of the atom to be assigned
     * @param res  The residue number to which the atom is assigned
     */
    void setResNum(int atom, int res);
    /**
     * Method to parse into a system
     */
    gcore::System parse();
    /**
     * Method to parse into an already existing system
     */
    void parse(gcore::System &sys);
    /**
     * Method to calculate the 1,4 neighbours from the bonds
     */
    void get14s();
    /**
     * Method that removes all references to all atoms with negative iac
     * it returns a vector of integers, with the new atom numbers, such that
     * these can also be used elsewhere
     */
    std::vector<int> removeAtoms();

    /**
     * Method that moves atoms forward or backward in the topology sequence
     * @param fromto list of pairs of old and new atom index which will be 
     *        applied sequentially
     */
    void moveAtoms(std::vector<std::pair<int, int> > fromto);

    // Accessors
    /**
     * Accessor, returns vector of atoms
     */
    std::vector<AtomTopology> & atoms();
    /**
     * Accessor, returns the vector of bonds
     */
    std::set<Bond> & bonds();
    /**
     * Accessor, returns the vector of constraints
     */
    std::set<Constraint> & constraints();
    /**
     * Accessor, returns the vector of dipole bonds
     */
    std::set<Bond> & dipoleBonds();
    /**
     * Accessor, returns the vector of angles
     */
    std::set<Angle> & angles();
    /**
     * Accessor, returns the vector of dihedrals
     */
    std::set<Dihedral> & dihedrals();
    /**
     * Accessor, returns the set of crossdihedrals
     */
    std::set<CrossDihedral> & crossdihedrals();
    /**
     * Accessor, returns the vector of impropers
     */
    std::set<Improper> & impropers();
    /**
     * Accessor, returns the vector of LJ exceptions
     */
    std::set<LJException> & ljexceptions();
    /**
     * Accessor, returns the vector with residue names
     */
    std::vector<std::string> & resNames();
    /**
     * Accessor, returns the map with residue numbers
     */
    std::map<int, int> & resMap();
  private:
    /**
     * Method that reduces the atoms array
     */
    void _reduceAtoms(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the residues
     */
    void _reduceResidues(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the bond vector
     */
    void _reduceBonds(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the constraint vector
     */
    void _reduceConstraints(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the dipole bond vector
     */
    void _reduceDipoleBonds(std::set<int> &rem, std::vector<int> &ren);    
    /**
     * Method that reduces the angle vector
     */
    void _reduceAngles(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the improper vector
     */
    void _reduceImpropers(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the dihedrals vector
     */
    void _reduceDihedrals(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the dihedrals vector
     */
    void _reduceCrossDihedrals(std::set<int> &rem, std::vector<int> &ren);
    /**
     * Method that reduces the LJ exception vector
     */
    void _reduceLJExceptions(std::set<int> &rem, std::vector<int> &ren);

  }; /* class LinearTopology */

} /* Namespace */

inline std::vector<gcore::AtomTopology> & gcore::LinearTopology::atoms() {
  return d_atom;
}

inline std::set<gcore::Bond> & gcore::LinearTopology::bonds() {
  return d_bond;
}

inline std::set<gcore::Constraint> & gcore::LinearTopology::constraints() {
  return d_constraints;
}

inline std::set<gcore::Bond> & gcore::LinearTopology::dipoleBonds() {
  return d_dipole_bond;
}

inline std::set<gcore::Angle> & gcore::LinearTopology::angles() {
  return d_angle;
}

inline std::set<gcore::Dihedral> & gcore::LinearTopology::dihedrals() {
  return d_dihedral;
}

inline std::set<gcore::CrossDihedral> & gcore::LinearTopology::crossdihedrals() {
  return d_crossdihedral;
}

inline std::set<gcore::Improper> & gcore::LinearTopology::impropers() {
  return d_improper;
}

inline std::set<gcore::LJException> & gcore::LinearTopology::ljexceptions() {
  return d_ljexception;
}

inline std::vector<std::string> & gcore::LinearTopology::resNames() {
  return d_resname;
}

inline std::map<int, int> & gcore::LinearTopology::resMap() {
  return d_resmap;
}
inline void gcore::LinearTopology::addAtom(const gcore::AtomTopology &a) {
  d_atom.push_back(a);
}

inline void gcore::LinearTopology::addBond(const gcore::Bond &b) {
  d_bond.insert(b);
}

inline void gcore::LinearTopology::addConstraint(const gcore::Constraint &b) {
  d_constraints.insert(b);
}

inline void gcore::LinearTopology::addDipoleBond(const gcore::Bond &b) {
  d_dipole_bond.insert(b);
}

inline void gcore::LinearTopology::addAngle(const gcore::Angle &a) {
  d_angle.insert(a);
}

inline void gcore::LinearTopology::addDihedral(const gcore::Dihedral &d) {
  d_dihedral.insert(d);
}

inline void gcore::LinearTopology::addCrossDihedral(const gcore::CrossDihedral &d) {
  d_crossdihedral.insert(d);
}

inline void gcore::LinearTopology::addImproper(const gcore::Improper &i) {
  d_improper.insert(i);
}

inline void gcore::LinearTopology::addLJException(const gcore::LJException &lj) {
  d_ljexception.insert(lj);
}

inline void gcore::LinearTopology::setResName(int i, const std::string &s) {
  if (int(d_resname.size()) <= i) d_resname.resize(i + 1);
  d_resname[i] = s;
}

inline void gcore::LinearTopology::setResNum(int i, int j) {
  d_resmap[i] = j;
}

#endif



