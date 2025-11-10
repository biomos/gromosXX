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

// fit_PositionUtils.h

#ifndef INCLUDED_FIT_PositionUtils
#define INCLUDED_FIT_PositionUtils


namespace utils
{
  class AtomSpecifier;
}

namespace gmath{
  class Vec;
  class Matrix;
}

namespace gcore{
  class System;
  class Molecule;
}

namespace fit{
  class Reference;
  /**
   * Class PositionUtils
   * Class that contains several operations to perform on your system
   *
   * Different functions that change the coordinates of your system are 
   * defined. It only works for the solute!
   *
   * @class PositionUtils
   * @author R. Buergi
   * @ingroup fit
   * @sa Reference
   */
  class PositionUtils{
    // not implemented
    PositionUtils (const PositionUtils&);
    PositionUtils &operator=(const PositionUtils &);
  public:
    /**
     * PositionUtils constructor
     */
    PositionUtils(){}
    /**
     * PositionUtils deconstructor
     */
    ~PositionUtils(){}
    
    /**
     * Method to calculat the centre of mass of your system
     * @return returns a Vec with the centre of mass position
     */
    static gmath::Vec com(const gcore::System &);
    /**
     * Method to calculate the centre of mass of your system, where a 
     * weight is added from the Reference (might be zero for many atoms)
     * @return a vector with the centre of mass position
     */
    static gmath::Vec com(const gcore::System &, const Reference &);
    /**
     * Method to calculate the centre of mass of your system, where a 
     * only considering atoms in atoms
     * @return a vector with the centre of mass position
     */
    static gmath::Vec com(const gcore::System &, utils::AtomSpecifier & atoms);
    /**
     * Method to calculate the velocity of the centre of mass of your system, where a 
     * only considering atoms in atoms
     * @return a vector with the velocity of the centre of mass
     */
    static gmath::Vec com_v(const gcore::System &, utils::AtomSpecifier & atoms);
    /**
     * Method to calculate the centre of geometry of your system
     * @return a Vec with the centre of geometry position
     */
    static gmath::Vec cog(const gcore::System &);
    /**
     * Method to calculate the centre of geometry of your system, where
     * a weight is added from the Reference (might be zero for many atoms)
     * @return a Vec with the centre of geometry position
     */
    static gmath::Vec cog(const gcore::System &, const Reference &);
    /**
     * Method to calculate the centre of geometry of your system, where
     * only atoms specified by atoms are considered
     * @return a Vec with the centre of geometry position
     */
    static gmath::Vec cog(const gcore::System &, utils::AtomSpecifier & atoms);
  
    /**
     * Method to move all solute atoms of your system by a Vec
     */
    static void translate(gcore::System *, const gmath::Vec &);
    /**
     * Method to move a solute molecule of your system by a Vec
     */
    static void translate(gcore::Molecule &, const gmath::Vec &);
    /**
     * Method to rotate all solute atoms of your system according to a 
     * rotation Matrix
     */
    static void rotate(gcore::System *, const gmath::Matrix &);
    /**
     * Method to rotate a solute molecule of your system according to a 
     * rotation Matrix
     */
    static void rotate(gcore::Molecule &, const gmath::Matrix &);
    /**
     * Method to calculate the matrix that rotates around 
     * the specified vector with the specified angle
     */
    static gmath::Matrix rotateAround(gmath::Vec v, double a);
    
    /**
     * Method to translate the System in such a way that its centre of
     * mass is at the origin
     */
    static gmath::Vec shiftToCom(gcore::System *);
    /** 
     * Method to translate the System in such a way that its centre of 
     * mass (calculated with the weights from the Reference) is at the 
     * origin
     */
    static gmath::Vec shiftToCom(gcore::System *, const Reference &);
    /**
     * Method to translate the System in such a way that its centre of 
     * mass (calculated using the atoms in the AtomSpecifier atoms) is at 
     * the origin
     */
    static gmath::Vec shiftToCom(gcore::System *, utils::AtomSpecifier & atoms);
    /**
     * Method to translate the System in such a way that its centre of
     * geometry is at the origin
     */
    static gmath::Vec shiftToCog(gcore::System *);
    /**
     * Method to translate the System in such a way that its centre of 
     * geometry (calculated with the weights from the Reference) is at 
     * the origin
     */
    static gmath::Vec shiftToCog(gcore::System *, const Reference &);
    /**
     * Method to translate the System in such a way that its centre of 
     * geometry (calculated using the atoms in the AtomSpecifier atoms) is at 
     * the origin
     */
    static gmath::Vec shiftToCog(gcore::System *, utils::AtomSpecifier & atoms);

    /**
     * Method to get maximum coordinates of the system 
     */
    static gmath::Vec getmaxcoordinates(gcore::System *, bool heavy);
    /**
     * Method to get minimum coordinate of the system 
     */
    static gmath::Vec getmincoordinates(gcore::System *, bool heavy);
  };

}

#endif
