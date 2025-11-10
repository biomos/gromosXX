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

// gcore_GromosForceField.h

#ifndef INCLUDED_GROMOSFORCEFIELD
#define INCLUDED_GROMOSFORCEFIELD

#include "LJExceptionType.h"
#include "AtomPair.h"
#include <map>
#include <string>

namespace gcore{

class GromosForceField_i;
class MassType;
class VirtualAtomType;
class BondType;
class AngleType;
class DihedralType;
class ImproperType;
class LJExceptionType;
class LJType;
class CGType;
/**
 * Class GromosForceField
 * Purpose: contains all force field parameters (ifp-file in gromos96)
 *
 * Description:
 * The GromosForceField contains all force field parameters that can be 
 * read from a topology. This roughly corresponds to the ifp-file in 
 * gromos96. Exceptions are the values of Fpepsi 
 * (@f$\frac{1}{4\pi\epsilon_0}@f$), Hbar (@f$\hbar@f$) and the speed of light
 * (spdl), which find a place in both the GromosForceField and the BuildingBlock
 * classes. This is due to the gromos96 structure
 *
 * @class GromosForceField
 * @author R. Buergi
 * @ingroup gcore
 * @sa gcore::MassType
 * @sa gcore::VirtualAtomType
 * @sa gcore::BondType
 * @sa gcore::AngleType
 * @sa gcore::DihedralType
 * @sa gcore::ImproperType
 * @sa gcore::LJType
 * @sa gcore::BuildingBlock
 */
class GromosForceField{
  GromosForceField_i *d_this;
  // not implemented
  GromosForceField &operator=(const GromosForceField &);

 public:
  /**
   * GromosForceField constructor
   */
  GromosForceField();
  /**
   * GromosForceField copy constructor
   */
  GromosForceField(const GromosForceField &);
  /**
   * GromosForceField deconstructor
   */
  ~GromosForceField();
  // Methods
  /**
   * Method to set the value of Fpepsi  (@f$\frac{1}{4\pi\epsilon_0}@f$)
   */
  void setFpepsi(double fpepsi);
  /**
   * Method to set the value of hbar (@f$\hbar@f$)
   */
  void setHbar(double hbar);
  /**
   * Method to set the value of the speed of light (spdl)
   */
  void setSpdl(double spdl);
  /**
   * Method to set the value of kB (@f$k_B@f$)
   */
  void setBoltz(double boltz);
  /**
   * Method to set the force field code
   */
  void setForceField(std::string code);
  /**
   * Method to add an Atom Type
   */
  void addAtomTypeName(const std::string &str);
  /**
   * Method to add a Mass Type
   * @param b MassType to add
   */
  void addMassType(const MassType &b);
  /**
   * Method to add a Virtual Atom Type
   * @param va VirtualAtomType to add
   */
  void addVirtualAtomType(const VirtualAtomType &va);
  /**
   * Method to add a Bond Type
   * @param b BondType to add
   */
  void addBondType(const BondType &b);
  /**
   * Method to add an Angle Type
   * @param b AngleType to add
   */
  void addAngleType(const AngleType &b);
  /**
   * Method to add a Dihedral Type
   * @param b DihedralType to add
   */
  void addDihedralType(const DihedralType &b);
  /**
   * Method to add an Improper Type
   * @param b ImproperType to add
   */
  void addImproperType(const ImproperType &b);
  /**
   * Method to add a LJ exception type
   * @param b LJ exception to add
   */
  void addLJExceptionType(const LJExceptionType &b);
  /**
   * Method to set a Lennard Jones exception interaction for a specific atom pair
   * @param p An AtomPair defined by their atom numbers of the topology
   * @param l The corresponding LJType containing the VDW parameters for 
   *          this AtomPair
   */
  void setLJException(const AtomPair &p, const LJExceptionType &l);
  /**
   * Method to set a Lennard Jones interaction for a specific atom pair
   * @param p An AtomPair defined by their Integer Atom Codes (iac's)
   * @param l The corresponding LJType containing the VDW parameters for 
   *          this AtomPair
   */
  void setLJType(const AtomPair &p, const LJType &l);
  /**
   * Method to set a coarse grain Lennard Jones interaction for a specific 
   * atom pair
   * @param p An AtomPair defined by their Integer Atom Codes (iac's)
   * @param l The corresponding CGType containing the coarse grain VDW 
   *          parameters for this AtomPair
   */
  void setCGType(const AtomPair &p, const CGType &l);

  // Accessors
  /**
   * Accessor, returns the value of Fpepsi 
   *           ( = @f$\frac{1}{4\pi\epsilon_0}@f$)
   */
  double fpepsi()const;
  /**
   * Accessor, returns the value of Hbar ( = @f$\hbar@f$)
   */
  double hbar()const;
  /**
   * Accessor, returns the value of the speed of light ( = spdl)
   */
  double spdl()const;
  /**
   * Accessor, returns the value of kB ( = @f$k_B@f$)
   */
  double boltz()const;
  /**
   * Accessor, returns the force field code
   */
  std::string ForceField()const;
  /**
   * Accessor, returns the number of Atom Type Names
   */
  int numAtomTypeNames()const;
  /**
   * Accessor, returns the name of the i-th atom type
   */
  const std::string &atomTypeName(int i) const;
  /** 
   * Accessor, returns the number of MassTypes
   */
  int numMassTypes()const;
  /**
   * Accessor, returns the i-th MassType
   */
  const MassType &massType(int i)const;
  /**
   * Method, returns the Mass from a MassType
   * @param i The gromos96 MassType
   * @return The corresponding Mass
   */
  double findMass(int i)const;
  /**
   * Method, returns the MassType from a Mass
   * @param mass the mass
   * @return The corresponding mass type or -1 if non was found.
   */
  int findMassType(double mass)const;
  /** 
   * Accessor, returns the number of VirtualAtomTypes
   */
  int numVirtualAtomTypes()const;
  /**
   * Accessor, returns the VirtualAtomType by type
   */
  const VirtualAtomType &virtualAtomType(int i) const;
  /**
   * Accessor, returns the i-th VirtualAtomType
   */
  const VirtualAtomType &virtualAtomTypeLine(int i) const;
  /**
   * Accessor, returns if the virtual atom type is already present
   */
  const bool findVirtualAtomType(const int i) const;
  /** 
   * Accessor, returns the number of BondTypes
   */
  int numBondTypes()const;
  /**
   * Accessor, returns the i-th BondType
   */
  const BondType &bondType(int i) const;
  /**
   * Accessor, returns the number of AngleTypes
   */
  int numAngleTypes()const;
  /**
   * Accessor, returns the i-th AngleType
   */
  const AngleType &angleType(int i) const;
  /**
   * Accessor, returns the number of DihedralTypes
   */
  int numDihedralTypes()const;
  /** 
   * Accessor, returns the i-th DihedralType
   */
  const DihedralType &dihedralType(int i) const;
  /**
   * Accessor, returns the number of ImproperTypes
   */
  int numImproperTypes()const;
  /**
   * Accessor, returns the i-th ImproperType
   */
  const ImproperType &improperType(int i) const;
  /**
   * Accessor, returns the number of LJ exception types
   */
  int numLJExceptionTypes()const;
  /**
   * Accessor, returns the i-th LJ Exception
   */
  const LJExceptionType &ljExceptionType(int i) const;
  /**
   * Accessor, returns the LJExceptionType for the specified AtomPair (gromos numbers!)
   */
  const std::map<AtomPair, LJExceptionType> &ljException() const;
  /**
   * Accessor, returns the number of LJTypes
   */
  int numLJTypes()const;
  /**
   * Accessor, returns the LJType for the specified AtomPair
   */
  const LJType &ljType(const AtomPair &p) const;
  /**
   * Accessor, returns the number of CGTypes
   */
  int numCGTypes()const;
  /**
   * Accessor, returns the cgType for the specified AtomPair
   */
  const CGType &cgType(const AtomPair &p) const;
  /**
   * Accessor, returns the dummy atom type or -1 if it is not found
   */
  int dummyAtomType()const;
  
};

}
#endif
