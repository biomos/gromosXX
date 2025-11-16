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

// utils_SimplePairlist.h
#ifndef INCLUDED_UTILS_SIMPLEPAIRLIST
#define INCLUDED_UTILS_SIMPLEPAIRLIST

#include "AtomSpecifier.h"

namespace bound{
  class Boundary;
}

namespace utils{
  class AtomSpecifier;
  /**
   * Class SimplePairlist
   *
   * A SimplePairlist is an AtomSpecifier, with additional functionality
   * to add all atoms within a certain cutoff from a reference atom.
   * Different cutoff schemes can be applied. And excluded atoms or 
   * 1,4 neighbours can be removed from the AtomSpecifier.
   * @class SimplePairlist
   * @author B.C. Oostenbrink
   * @ingroup utils
   */
  class SimplePairlist: public AtomSpecifier{
    /**
     * A pointer to the Periodic Boundary Conditions. Distances that
     * between atoms are always calculated over the nearestImage.
     */
    bound::Boundary *d_pbc;
    /**
     * The square of the cutoff
     */
    double d_cut2;
    /** 
     * A boolean which determines the Pairlist scheme that is employed
     */
    bool d_chargeGroupBased;
    /**
     * The reference atom
     */
    SpecAtom *d_atom;
    
  public:
    /**
     * Constructor
     */
    SimplePairlist(){};
    /**
     * Constructor
     *
     * @param sys The system for which the SimplePairlist
     *                          will be constructed
     * @param pbc The Periodic Boundary Conditions that need
     *            to be applied to the system
     * @param c The cutoff.
     */
    SimplePairlist(gcore::System &sys, bound::Boundary &pbc, double c);
    /**
     * A function to set the Periodic Boundary Conditions
     */
    void setPbc(bound::Boundary &pbc);
    /**
     * A function to set the cutoff
     */
    void setCutOff(double c);
    /**
     * A function to set the pairlist scheme
     * @param s Can take two values: "ATOMIC" will use an atomic 
     *        cutoff scheme, "CHARGEGROUP" will use a charge 
     *        group based cutoff scheme.
     */
    void setType(std::string s);
    /**
     * A function to set the reference atom for which the SimplePairlist is 
     * calculated
     * @param m the molecule number
     * @param a the atom number
     */
    void setAtom(int m, int a);
    /**
     * A function to set the reference atom for which the SimplePairlist is
     * calculated
     * @param  s an AtomSpecifier containing 1 atom
     */
    void setAtom(SpecAtom &s);
    /**
     * Calculates the SimplePairlist according to the scheme which has been 
     * set by setType();
     */
    void calc();
    /**
     * Calculates a charge group based atomic pairlist. Excluded atoms and 
     * 1,4 neighbours will be included. The system is first
     * gathered so that all charge groups are connected.
     */
    void calcCgb();
    /**
     * Calculates an atom based atomic pairlist. Excluded atoms and 1,4 
     * neighbours will be included.
     */
    void calcAtomic();
    /**
     * Calculates the SimplePairlist according to the scheme which has been 
     * set by setType();
     */
    void calc(const AtomSpecifier &B, double cutmin = 0.0);
    /**
     * Calculates a charge group based atomic pairlist. Excluded atoms and 
     * 1,4 neighbours will be included. The system is first
     * gathered so that all charge groups are connected.
     */
    void calcCgb(const AtomSpecifier &B, double cutmin = 0.0);
    /**
     * Calculates an atom based atomic pairlist. Excluded atoms and 1,4 
     * neighbours will be included. Only atoms within the specifiers B are
     * considered. The function assums B contains complete charge groups!
     */
    void calcAtomic(const AtomSpecifier &B, double cutmin = 0.0);
    /**
     * Function to remove all excluded atoms from the the SimplePairlist
     */
    void removeExclusions();
    /**
     * Function to remove all 1,4 neighbours from the SimplePairlist
     */
    void remove14Exclusions();
  protected:
    /**
     * Function to calculate the position of a charge group to which the 
     * specified atom belongs. For solute this is the centre of geometry of 
     * all atoms beloning to the charge group, for solvent it is the position 
     * of the first atom of the solvent molecule.
     * @param m the molecule number
     * @param a the atom number
     */
    gmath::Vec chargeGroupPosition(int m, int a);
    /**
     * Function to calculate the position of a charge group to which the
     * specified atom belongs. For solute this is the centre of geometry of all
     * atoms belonging to the charge group, for solvent it is the position of
     * the first atom of the solvent molecule. For a virtual atom, it is the 
     * position itself
     * @param s the atom
     */
    gmath::Vec chargeGroupPosition(SpecAtom &s);
  };
}

#endif 
