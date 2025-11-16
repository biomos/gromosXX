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

// gcore_SolventTopology.h

#ifndef INCLUDED_GCORE_SOLVENTTOPOLOGY
#define INCLUDED_GCORE_SOLVENTTOPOLOGY

#include <string>

namespace gcore{

  class SolventTopology_i;
  class GromosForceField;
  class AtomTopology;
  class Constraint;
  class ConstraintIterator;
  /**
   * Class SolventTopology
   * Purpose: Contains all topological information for a solvent Molecule
   *
   * Description:
   * The SolventTopology contains all topological information for a solvent 
   * molecule. In gromos a solvent molecule is treates different from a 
   * regular molecule: All atoms are excluded and in stead of bonds, angles
   * etc. the molecule is kept rigid by a set of distance restraints.<br>
   * For the Atoms there is a direct accessor to the AtomTopologies, but 
   * these Constraints you access via a ConstraintIterator.
   *
   * @class SolventTopology
   * @author M. Kastenholz and C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::Solvent
   * @sa gcore::Constraint
   */
  class SolventTopology{
    SolventTopology_i *d_this;
    // This class contains all topological information
    /**
     * Constraint Iterator for the a SolventTopology
     *
     * The SolventTopology Constraint iterator is used to loop over the  
     * Constraints in a SolventTopology. 
     * It is constructed with the SolventTopology as an argument. Use the 
     * ++ operator to move to the next Constraint. The () operator returns 
     * the current Constraint. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Constraint list.
     * @author M. Kastenholz and C. Oostenbrink
     */

    friend class ConstraintIterator;
    
    // not implemented
    SolventTopology &operator=(const SolventTopology &);

  public:
    /**
     * SolventTopology constructor
     */
    SolventTopology();
    /**
     * SolventTopology copy constructor
     * @param & SolventTopology to be copied
     */
    SolventTopology(const SolventTopology &);
    /**
     * SolventTopology deconstructor
     */
    ~SolventTopology();
    
    // Methods
    /**
     * Method to add an atom to the SolventTopology
     * @param a An AtomTopology that is to be added; note that even though
     *          in gromos the chargeGroup and Exclusions are not used in 
     *          solvent, you can still have these in your AtomTopology
     */
    void addAtom(const AtomTopology &a);
    /**
     * Method to add a constraint to the SolventTopology
     * @param b The Constraint that is to be added
     */
    void addConstraint(const Constraint &b);
    /**
     * Method to set the name of the solvent
     */
    void setSolvName(const std::string &s);
    /**
     * Method to determine which atoms are hydrogens based on the mass
     */
    void setHmass(double mass);
    /**
     * Method to determine which atoms are hydrogens based on the iac
     */
    void setHiac(int iac);
    /**
     * Method to clear all isH flags of the atoms
     */
    void clearH();
    
    // Accessors
    /**
     * Accessor, returns the number of atoms in the SolventTopology. Note 
     * that this is the number of atoms in a single solvent molecule.
     */
    int numAtoms()const;
    /**
     * Accessor, returns an AtomTopology for the i-th atom in a solvent 
     * molecule as a const
     */
    const AtomTopology& atom(int i) const; 
    /**
     * Accessor, returns an AtomTopology for the i-th atom in a solvent
     * molecule as a const
     */
    AtomTopology & atom(int i);
    /**
     * Accessor, returns the name of the solvent
     */
    const std::string &solvName()const;
    /**
     * the number of constraints
     */
    int numConstraints()const;
    
  }; /* class MoleculeTopology */


  class ConstraintIterator_i;

  class ConstraintIterator{
    ConstraintIterator_i *d_this;
    // not implemented
    ConstraintIterator();
    ConstraintIterator(const ConstraintIterator&);
    ConstraintIterator &operator=(const ConstraintIterator &);
  public:
    ConstraintIterator(const SolventTopology &mt);
    ~ConstraintIterator();
    void operator++();
    const Constraint &operator()()const;
    operator bool()const;
  };

} /* Namespace */ 
#endif



