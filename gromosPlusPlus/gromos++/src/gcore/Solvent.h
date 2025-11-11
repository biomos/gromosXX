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

// gcore_Solvent.h
/**
 * Class Solvent
 * Addition: velocity configuration added to Molecule definition;
 * Author  : gee          
 */

#ifndef INCLUDED_GCORE_SOLVENT
#define INCLUDED_GCORE_SOLVENT

#include <vector>
#include <cassert>

namespace gmath{
class Vec;
}

using gmath::Vec;

namespace gcore{

  class SolventTopology;

  /** 
   * Class Solvent
   * Purpose: Contains coordinates and topological information of the solvent
   *
   * Description:
   * The gromos++ Solvent is a special kind of molecule. The class contains
   * topological information for the solvent (one) as well as the coordinates
   * for ALL solvent molecules in the system (many molecules). This of
   * course because it is not very useful to store 15000 copies of the same 
   * topology.
   *
   * @class Solvent
   * @author M. Kastenholz and C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::SolventTopology
   * @sa gcore::System
   */
  class Solvent{
    SolventTopology *d_mt;
    std::vector<Vec*> d_pos;
    int d_numPos;
    std::vector<Vec*> d_vel;
    int d_numVel;
    std::vector<Vec*> d_cosDisplacement;
    int d_numCosDisplacements;

    // not implemented
    Solvent();
    //Solvent &operator=(const Solvent &);

  public:
    /**
     * Solvent Constructor
     * @param mt a SolventTopology, this creates a solvent without coordinates
     */
    Solvent(const SolventTopology &mt);
    /**
     * Solvent Copy Constructor
     * @param & Solvent to be copied
     */
    Solvent(const Solvent &);
    /**
     * Solvent deconstructor
     */
    ~Solvent();
    
    // Methods
    /**
     * Member operator = copies one Solvent into the other
     */
    Solvent &operator=(const Solvent &s);
     /**
     * Accessor, returns a SolventTopology containing the topological 
     * information for one (1) solvent molecule
     */
    SolventTopology &topology(); 
    /**
     * Accessor, returns a pointer to a vector with the coordinates of 
     * solvent atom i. 
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    Vec &pos(int i);
    /**
     * Accessor, returns a pointer to a vector with the coordinates of 
     * solvent atom i as a const
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    const Vec &pos(int i)const;
    /**
     * Accessor, returns a pointer to a vector with the velocity coordinates
     * of solvent atom i. 
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    Vec &vel(int i);
    /**
     * Accessor, returns a pointer to a vector with the velocity coordinates
     * of solvent atom i as a const
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    const Vec &vel(int i)const;
    /**
     * Accessor, returns a pointer to a vector with the charge-on-spring
     * displacement of solvent atom i. 
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    Vec &cosDisplacement(int i);
    /**
     * Accessor, returns a pointer to a vector with the charge-on-spring
     * displacement of solvent atom i as a const. 
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    const Vec &cosDisplacement(int i)const;
    /**
     * Accessor, returns a SolventTopology containing the topological 
     * information for one (1) solvent molecule
     */
    const SolventTopology &topology() const; 
    /**
     * Method to add an atom to the Solvent.
     *
     * The user is responsible for adding complete molecules, that is adding
     * all atoms of solvent molecule to the Solvent (and in the correct order)
     * @param v A gmath::Vec containing three coordinates
     */
    void addPos(Vec v);
    /**
     * Method to rescale the number of atoms in the Solvent
     *
     * This method allows you to set the total number of solvent atoms we have
     * If i < numCoords then the memory containing the coordinates for all 
     * atoms >i is released.
     * @param i The new total number of Solvent atoms in the class
     */
    void setNumPos(int i);
    /**
     * Accessor, returns the number of atoms for which coordinates are 
     * stored in the Solvent class. So this is the total number of solvent
     * atoms the class knows about (=number of atoms per solvent molecules
     * * number of solvent molecules)
     */
    int numPos()const;
    /**
     * Method to add a velocity vector to the Solvent velocity configuration.
     *
     * The user is responsible for adding complete molecules, that is adding
     * velocity vectors for all atoms of solvent molecule to the Solvent
     * (and in the correct order)
     * @param v A gmath::Vec containing three coordinates
     */
    void addVel(Vec v);
    /**
     * Method to rescale the number of velocity coordinates in the Solvent
     * velocity configuration.
     * This method allows you to set the total number of solvent atoms we have
     * If i < numCoords then the memory containing the coordinates for all 
     * atoms >i is released.
     * @param i The new total number of Solvent atoms in the class
     */
    void setNumVel(int i);
    /**
     * Accessor, returns the number of atoms for which velocity coordinates
     * are stored in the Solvent class. So this is the total number of solvent
     * atom velocity vectors the class knows about
     * (=number of atoms per solvent molecule * number of solvent molecules)
     */
    int numVel()const;
     /**
     * Method to add a COS displacement vector to the Solvent configuration.
     *
     * The user is responsible for adding complete molecules, that is adding
     * COS displacement vectors for all atoms of solvent molecule to the 
     * Solvent (and in the correct order)
     * @param c A gmath::Vec containing three coordinates
     */
    void addCosDisplacement(Vec c);
    /**
     * Method to rescale the number of COS displacements in the Solvent
     * configuration.
     * This method allows you to set the total number of solvent atoms we have
     * If i < numCoords then the memory containing the coordinates for all 
     * atoms >i is released.
     * @param i The new total number of Solvent atoms in the class
     */
    void setNumCosDisplacements(int i);
    /**
     * Accessor, returns the number of atoms for which COS displacements
     * are stored in the Solvent class. So this is the total number of solvent
     * atom COS displacement vectors the class knows about
     * (=number of atoms per solvent molecule * number of solvent molecules)
     */
    int numCosDisplacements()const;
    /**
     * Accessor, returns the number of atoms as the maximum of either
     * numPos(), numVel() or numCosDisplacements()
     */
    int numAtoms()const;
    
    
  }; /* class Solvent */

  inline Vec &Solvent::pos(int i){
    assert(i < (this->numPos()));
    return *d_pos[i];
  }

  inline const Vec &Solvent::pos(int i)const{
    assert (i < (this->numPos()));
    return *d_pos[i];
  }

  inline Vec &Solvent::vel(int i){
    assert(i < (this->numVel()));
    return *d_vel[i];
  }

  inline const Vec &Solvent::vel(int i)const{
    assert (i < (this->numVel()));
    return *d_vel[i];
  }
  
  inline Vec &Solvent::cosDisplacement(int i){
    assert(i < (this->numCosDisplacements()));
    return *d_cosDisplacement[i];
  }

  inline const Vec &Solvent::cosDisplacement(int i)const{
    assert (i < (this->numCosDisplacements()));
    return *d_cosDisplacement[i];
  }
  
} /* Namespace */ 
#endif

