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

// gcore_MoleculeTopology.h

#ifndef INCLUDED_GCORE_MOLECULETOPOLOGY
#define INCLUDED_GCORE_MOLECULETOPOLOGY

#include <set>
#include <string>

namespace gcore{

  class MoleculeTopology_i;
  class GromosForceField;
  class AtomTopology;
  class Bond;
  class Constraint;
  class Angle;
  class Dihedral;
  class CrossDihedral;
  class Improper;
  class LJException;
  class BondIterator;
  class BondDipoleIterator;
  class AngleIterator;
  class ImproperIterator;
  class DihedralIterator;
  class LJExceptionIterator;
  /**
   * Class MoleculeTopology
   * Purpose: Contains all topological information for a Molecule
   *
   * Description:
   * The MoleculeTopology contains all topological information for a 
   * Molecule. For the Atoms there is a direct accessor to the 
   * AtomTopologies. Bonds, Angles etc. you access via the BondIterator, 
   * AngleIterator, ImproperIterator, DihedralIterator and CrossDihedralIterator.
   *
   * @class MoleculeTopology
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::Molecule
   * @sa gcore::AtomTopology
   * @sa gcore::Bond
   * @sa gcore::Angle
   * @sa gcore::Improper
   * @sa gcore::Dihedral
   * @sa gcore::CrossDihedral
   */
  class MoleculeTopology{
    MoleculeTopology_i *d_this;
    // This class contains all topological information
    /**
     * Atom Iterator for a MoleculeTopology
     *
     * The MoleculeTopology Atom iterator is used to loop over the AtomTopologies 
     * in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next atom topology. The () operator returns the 
     * current atom topology. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the atom list.
     * @author R. Buergi
     */
    friend class AtomIterator;
    /**
     * Bond Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology Bond iterator is used to loop over the Bonds 
     * in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next Bond. The () operator returns the 
     * current Bond. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the bond list.
     * @author R. Buergi
     */
    friend class BondIterator;
    /**
     * Dipole Bond Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology BondDipole iterator is used to loop over the Dipole 
     * Bonds in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next Bond. The () operator returns the 
     * current Bond. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the bond list.
     * @author Victor H. Rusu
     */
    friend class BondDipoleIterator;    
    /**
     * Angle Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology Angle iterator is used to loop over the Angles 
     * in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next Angle. The () operator returns the 
     * current Angle. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the angle list.
     * @author R. Buergi
     */
    friend class AngleIterator;
    /**
     * Improper Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology Improper iterator is used to loop over the 
     * Impropers in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next Improper. The () operator returns the 
     * current Improper. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the improper list.
     * @author R. Buergi
     */
    friend class ImproperIterator;
    /**
     * Dihedral Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology Dihedral iterator is used to loop over the 
     * Dihedrals in a MoleculeTopology. 
     * It is constructed with the MoleculeTopology as an argument. Use the 
     * ++ operator to move to the next Dihedral. The () operator returns the 
     * current Dihedral. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Dihedral list.
     * @author R. Buergi
     */
    friend class DihedralIterator;

    /**
     * CrossDihedral Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology CrossDihedral iterator is used to loop over the
     * CrossDihedrals in a MoleculeTopology.
     * It is constructed with the MoleculeTopology as an argument. Use the
     * ++ operator to move to the next CrossDihedral. The () operator returns the
     * current CrossDihedral.
     * This can also be used as a boolean: the bool() returns 1 as long as
     * the iterator is not at the end of the CrossDihedral list.
     * @author N. Schmid
     */
    friend class CrossDihedralIterator;

    /**
     * LJException Iterator for the a MoleculeTopology
     *
     * The MoleculeTopology LJException iterator is used to loop over the
     * LJ Exceptions in a MoleculeTopology.
     * It is constructed with the MoleculeTopology as an argument. Use the
     * ++ operator to move to the next LJ exception. The () operator returns the
     * current LJ exception.
     * This can also be used as a boolean: the bool() returns 1 as long as
     * the iterator is not at the end of the LJException list.
     * @author A.P. Eichenberger
     */
    friend class LJExceptionIterator;
    


  public:
    /**
     * MoleculeTopology constructor
     */
    MoleculeTopology();
    /**
     * MoleculeTopology copy constructor
     * @param & MoleculeTopology to be copied
     */
    MoleculeTopology(const MoleculeTopology &);
    /**
     * MoleculeTopology deconstructor
     */
    ~MoleculeTopology();
    
    /**
     * MoleculeTopology = operator
     */
    MoleculeTopology &operator=(const MoleculeTopology &); 

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
     * @param b The Dipole Bond that is to be added; should be complete already
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
     * Accessor, returns the number of atoms in the MoleculeTopology
     */
    int numAtoms()const;
    /**
     * Accessor, returns the number of bonds in the MoleculeTopology
     */
    int numBonds()const;
    /**
     * Accessor, returns the number of constraints in the MoleculeTopology
     */
    int numConstraints()const;
    /**
     * Accessor, returns the number of dipole bonds in the MoleculeTopology
     */
    int numDipoleBonds()const;
    /**
     * Accessor, returns the number of angles in the MoleculeTopology
     */
    int numAngles()const;
    /**
     * Accessor, returns the number of impropers in the MoleculeTopology
     */
    int numImpropers()const;
    /**
     * Accessor, returns the number of dihedrals in the MoleculeTopology
     */
    int numDihedrals()const;
    /**
     * Accessor, returns the number of cross-dihedrals in the MoleculeTopology
     */
    int numCrossDihedrals()const;
    /**
     * Accessor, returns the number of LJ exception in the MoleculeTopology
     */
    int numLJExceptions()const;
    /**
     * Accessor, return the AtomTopology of the i-th atom in the 
     * MoleculeTopology
     */
    AtomTopology & atom(int i);
    /**
     * Accessor, return the AtomTopology of the i-th atom in the 
     * MoleculeTopology as a const
     */
    const AtomTopology& atom(int i) const; 
    /**
     * Accessor, returns the number of residues in the molecule
     */
    int numRes()const;
    /**
     * Accessor, returns the residue number to which the atom belongs
     * @param atom the atom number of an atom in the topology
     * @return the number of the residue to which this atom belongs
     */
    int resNum(int atom) const;
    /**
     * Accessor, returns the name of a residue number
     * @param i The number of a residue in the MoleculeTopology
     * @return A string with the name of the residue
     */
    const std::string &resName(int i)const;

    /**
     * Accessor to the set of constraints
     */
     std::set<Constraint> & constraints()const;


  }; /* class MoleculeTopology */

  class AtomIterator_i;
  class BondIterator_i;
  class BondDipoleIterator_i;
  class AngleIterator_i;
  class ImproperIterator_i;
  class DihedralIterator_i;
  class CrossDihedralIterator_i;
  class LJExceptionIterator_i;


  class AtomIterator{
    AtomIterator_i *d_this;
    // not implemented
    AtomIterator();
    AtomIterator(const AtomIterator&);
    AtomIterator &operator=(const AtomIterator &);
  public:
    AtomIterator(const MoleculeTopology &mt);
    ~AtomIterator();
    void operator++();
    const AtomTopology &operator()()const;
    AtomTopology & operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

  class BondIterator{
    BondIterator_i *d_this;
    // not implemented
    BondIterator();
    BondIterator(const BondIterator&);
    BondIterator &operator=(const BondIterator &);
  public:
    BondIterator(const MoleculeTopology &mt);
    ~BondIterator();
    void operator++();
    const Bond &operator()()const;
    Bond & operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };
  
  class BondDipoleIterator{
    BondDipoleIterator_i *d_this;
    // not implemented
    BondDipoleIterator();
    BondDipoleIterator(const BondDipoleIterator&);
    BondDipoleIterator &operator=(const BondDipoleIterator &);
  public:
    BondDipoleIterator(const MoleculeTopology &mt);
    ~BondDipoleIterator();
    void operator++();
    const Bond &operator()()const;
    Bond & operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };  

  class AngleIterator{
    AngleIterator_i *d_this;
    // not implemented
    AngleIterator();
    AngleIterator(const AngleIterator&);
    AngleIterator &operator=(const AngleIterator &);
  public:
    AngleIterator(const MoleculeTopology &mt);
    ~AngleIterator();
    void operator++();
    const Angle &operator()()const;
    Angle &operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

  class ImproperIterator{
    ImproperIterator_i *d_this;
    // not implemented
    ImproperIterator();
    ImproperIterator(const ImproperIterator&);
    ImproperIterator &operator=(const ImproperIterator &);
  public:
    ImproperIterator(const MoleculeTopology &mt);
    ~ImproperIterator();
    void operator++();
    const Improper &operator()()const;
    Improper &operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

  class DihedralIterator{
    DihedralIterator_i *d_this;
    // not implemented
    DihedralIterator();
    DihedralIterator(const DihedralIterator&);
    DihedralIterator &operator=(const DihedralIterator &);
  public:
    DihedralIterator(const MoleculeTopology &mt);
    ~DihedralIterator();
    void operator++();
    const Dihedral &operator()()const;
    Dihedral &operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

  class CrossDihedralIterator{
    CrossDihedralIterator_i *d_this;
    // not implemented
    CrossDihedralIterator();
    CrossDihedralIterator(const CrossDihedralIterator&);
    CrossDihedralIterator &operator=(const CrossDihedralIterator &);
  public:
    CrossDihedralIterator(const MoleculeTopology &mt);
    ~CrossDihedralIterator();
    void operator++();
    const CrossDihedral &operator()()const;
    CrossDihedral &operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

  class LJExceptionIterator{
    LJExceptionIterator_i *d_this;
    // not implemented
    LJExceptionIterator();
    LJExceptionIterator(const LJExceptionIterator&);
    LJExceptionIterator &operator=(const LJExceptionIterator &);
  public:
    LJExceptionIterator(const MoleculeTopology &mt);
    ~LJExceptionIterator();
    void operator++();
    const LJException &operator()()const;
    LJException &operator()();
    operator bool()const;
    bool last()const;
    bool first()const;
  };

} /* Namespace */ 
#endif



