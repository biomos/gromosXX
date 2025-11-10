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

// gcore_System.h

#ifndef INCLUDED_GCORE_SYSTEM
#define INCLUDED_GCORE_SYSTEM

#include <vector>
#include <cassert>

#include "VirtualAtoms.h"

namespace gmath{
class Vec;
}
using gmath::Vec;

namespace gcore{

  class Solvent;
  class Molecule;
  class VirtualAtoms;
  class Box;
  class Remd;
  class Weight;
  /**
   * Class System
   * The system class in gromos++ contains everything for the 
   * that you were interested in. Coordinates and topological information 
   * of the Molecules and Solvent.
   *
   * Description:
   * The System class contains all information of the molecules you have
   * or will simulate. Coordinates (including box sizes) and topological 
   * information for your solute and solvent molecules are accessible 
   * from the system.
   *
   * @class System
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::Molecule
   * @sa gcore::Solvent
   * @sa gcore::Box
   * @sa gcore::Remd
   */
  class System{
    std::vector<Molecule*> d_mol;
    std::vector<Solvent*> d_sol;
    VirtualAtoms d_vas;
    std::vector<int*> d_temperatureGroup;
    std::vector<int*> d_pressureGroup;
    Box *d_box;
    Remd *d_remd;
    Weight *d_weight;

  public:
    //Constructors
    /**
     * System Constructor
     */
    System();
    /**
     * System copy constructor
     */
    System(const System &sys);
    /**
     * System deconstructor
     */
    ~System();
    /**
     * Boolean to indicate whether a Box block has been read in.
     */
    bool hasBox;
    /**
     * Boolean to indicate whether a Coordinate block has been read in.
     */
    bool hasPos;
    /**
     * Boolean to indicate whether a lattice shift block has been read in.
     */
    bool hasLatticeshifts;
    /**
     * Boolean to indicate whether a Velocity block has been read in.
     */
    bool hasVel;
    /**
     * Boolean to indicate whether a COSDISPLACEMENT block has been read in.
     */
    bool hasCosDisplacements;
    /**
     * Boolean to indicate whether a REMD block has been read in.
     */
    bool hasRemd;
    
    // Methods
    /**
     * Member operator = copies one System into the other
     */
    System &operator=(const System &sys);
    /**
     * Method to add a Molecule to your system
     *
     * After adding a molecule to your system, you can still modify
     * the coordinates, but the topological information is now considered
     * as const.
     * @param mol Molecule to be added
     */
    void addMolecule(const Molecule &mol);
    /**
     * Method to add a Solvent to your system
     *
     * After adding a solvent to your system, you can still modify the 
     * coordinates, but the topological information is now considered as
     * const.<br>
     * Note that even though gromos96 does not support this option, gromos++
     * can in principle handle systems with multiple solvents.
     * @param sol Solvent to be added
     */
    void addSolvent(const Solvent &sol);
    /**
     * Method to add a temperature group to your system
     * @param tg Last atom number in temperature group to add
     */
    void addTemperatureGroup(const int &tg);
    /**
     * Method to set a pressure group to your system
     * @param pg Last atom number in pressure group to add
     */
    void addPressureGroup(const int &pg);
    /**
     * Method to add a set of virtual atoms
     */
    void addVirtualAtoms(gcore::VirtualAtoms &vas);
    /**
     * Method to add a single virtual atom
     */
    void addVirtualAtom(std::vector<int> conf, int type, double dish, double disc, int iac, double charge, gcore::Exclusion e, gcore::Exclusion e14);
 
    // Accessors
    /** 
     * Accessor, returns the i-th Molecule in the System as a const
     */
    const Molecule &mol(int i)const;
    /**
     * Accessor, returns the i-th Molecule in the System
     */
    Molecule &mol(int i);
    /**
     * Accessor, returns the i-th Solvent in the System as a const. Remember
     * that this is the i-th solvent type. All solvent molecules of this
     * type are stored within the Solvent class
     */
    const Solvent &sol(int i)const;
    /**
     * Accessor, returns the i-th Solvent in the System. Remember
     * that this is the i-th solvent type. All solvent molecules of this
     * type are stored within the Solvent class
     */
    Solvent &sol(int i);
    /**
     * Accessor, returns the Box dimensions of the System as a const
     */
    const Box &box()const;
    /** 
     * Accessor, returns the Box dimensions of the System
     */
    Box &box();
    /**
     * Accessor, returns the REMD information
     */
    Remd &remd();
    /**
     * Accessor, returns the REMD information as a const
     */
    Remd &remd()const;
    /**
     * Accessor, returns the last atom number of the i-th temperature group
     */
    int &temperatureGroup(int i);
    /**
     * Accessor, returns the last atom number of the i-th temperature group as a const
     */
    const int &temperatureGroup(int i) const;
    /**
     * Accessor, returns the last atom number of the i-th pressure group
     */
    int &pressureGroup(int i);
    /**
     * Accessor, returns the last atom number of the i-th pressure group as a const
     */
    const int &pressureGroup(int i) const;
    
    /**
     * Accessor, returns the VirtualAtoms of the System as a const
     */
    const VirtualAtoms vas()const;
    /** 
     * Accessor, returns the VirtualAtoms of the System
     */
    VirtualAtoms vas();
    /**
     * Accessor, returns the number of Molecules in the System
     */
    int numMolecules()const;
    /**
     * Accessor, returns the number of Solvents in the System. Note 
     * that this again the number of different Solvent types in the 
     * system. The number of solvent molecules (or rather total number
     * of solvent atoms) can be obtained from the Solvent class.
     */
    int numSolvents()const;
    /**
     * Accessor, returns the number of temperature groups in the System
     */
    int numTemperatureGroups()const;
    /**
     * Accessor, returns the number of temperature groups in the System
     */
    int numPressureGroups()const;
     /**
     * Accessor, returns a list of primary atoms in the System
     */
    mutable std::vector<Vec> primlist;
    /**
     * Accessor to the weight of the configuration
     */
    Weight & weight();
    /**
     * Accessor to the weight of the configuration (const version)
     */
    const Weight & weight() const;
  };

  inline const Molecule &System::mol(int i)const{
    assert (i < int(this->numMolecules()));
    return *d_mol[i];
  }

  inline Molecule &System::mol(int i){
    assert (i < int(this->numMolecules()));
    return *d_mol[i];
  }

  inline const Solvent &System::sol(int i)const{
    assert (i < this->numSolvents());
    return *d_sol[i];
  }

  inline Solvent &System::sol(int i){
    assert (i < this->numSolvents());
    return *d_sol[i];
  }

  inline const int &System::temperatureGroup(int i)const{
    assert (i < this->numTemperatureGroups());
    return *d_temperatureGroup[i];
  }

  inline int &System::temperatureGroup(int i){
    assert (i < this->numTemperatureGroups());
    return *d_temperatureGroup[i];
  }

  inline const int &System::pressureGroup(int i)const{
    assert (i < this->numPressureGroups());
    return *d_pressureGroup[i];
  }

  inline int &System::pressureGroup(int i){
    assert (i < this->numPressureGroups());
    return *d_pressureGroup[i];
  }

  inline const Box &System::box()const{
    return *d_box;
  }

  inline Box &System::box(){
    return *d_box;
  }
  inline Remd &System::remd(){
    return *d_remd;
  }
  inline Remd &System::remd()const
  {
    return *d_remd;
  }
  
  inline const VirtualAtoms System::vas()const{
    return d_vas;
  }

  inline VirtualAtoms System::vas(){
    return d_vas;
  }
  inline int System::numMolecules()const{
    return d_mol.size();
  }

  inline int System::numSolvents()const{
      return d_sol.size();
  }

  inline int System::numTemperatureGroups()const{
    return d_temperatureGroup.size();
  }

  inline int System::numPressureGroups()const{
      return d_pressureGroup.size();
  }
  
  inline Weight & System::weight() {
    return *d_weight;
  }
  
  inline const Weight & System::weight() const {
    return *d_weight;
  }
  
}
#endif
