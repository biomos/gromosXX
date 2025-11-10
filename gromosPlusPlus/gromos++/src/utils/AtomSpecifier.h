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

/**
 * @file AtomSpecifier.h
 * AtomSpecifier methods
 */

// Class that contains a sequential list of specific atoms

#ifndef INCLUDED_UTILS_ATOMSPECIFIER
#define INCLUDED_UTILS_ATOMSPECIFIER

#include <vector>
#include <string>

#include "../gromos/Exception.h"


// minimal complete headers
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"

#include "../utils/VirtualAtom.h"

namespace gmath
{
  class Vec;
}

namespace utils
{
  enum spec_type{ spec_solute, spec_solvent, spec_virtual };

  class AtomSpecifier;

  /**
   * Class SpecAtom
   * purpose: interface to access a specific atom (in an AtomSpecifier)
   * unified access to solute, solvent and virtual atoms.
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class SpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class SpecAtom
  {
  public:
    SpecAtom(){}

    SpecAtom(gcore::System &sys, int m, int a) : d_sys(&sys), d_mol(m), d_atom(a) {}
    SpecAtom(SpecAtom const & s) : d_sys(s.d_sys), d_mol(s.d_mol), d_atom(s.d_atom) {}

    virtual ~SpecAtom() {};

    virtual spec_type type() const { return spec_solute; }

    virtual SpecAtom * clone() { return new SpecAtom(*this); }

    virtual int mol() const { return d_mol; }
    virtual int atom()const { return d_atom; }

    virtual gmath::Vec & pos() { return d_sys->mol(d_mol).pos(d_atom); }
    virtual gmath::Vec const & pos() const { return d_sys->mol(d_mol).pos(d_atom); }
    virtual gmath::Vec & cosDisplacement() { return d_sys->mol(d_mol).cosDisplacement(d_atom); }
    virtual gmath::Vec const & cosDisplacement() const { return d_sys->mol(d_mol).cosDisplacement(d_atom); }
    virtual gmath::Vec & vel() { return d_sys->mol(d_mol).vel(d_atom); }
    virtual gmath::Vec const & vel() const { return d_sys->mol(d_mol).vel(d_atom); }
    virtual std::string name()const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).name();
    }

    virtual int iac() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).iac();
    }
    virtual double radius() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).radius();
    }
    virtual double charge() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).charge();
    }
    virtual bool isPolarisable() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).isPolarisable();
    }
    virtual double poloffsiteGamma() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).poloffsiteGamma();
    }
    virtual int poloffsiteI() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).poloffsiteI();
    }
    virtual int poloffsiteJ() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).poloffsiteJ();
    }
    virtual double cosCharge() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).cosCharge();
    }
    virtual double mass() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).mass();
    }
    virtual int resnum() const
    {
      return d_sys->mol(d_mol).topology().resNum(d_atom);
    }
    virtual std::string resname() const
    {
      return d_sys->mol(d_mol).topology().resName(resnum());
    }

    virtual void setSystem(gcore::System &sys)
    {
      d_sys = &sys;
    }

    virtual std::string toString() const;

    virtual utils::AtomSpecifier conf();

    virtual int virtualType() { return 0; }


  protected:
    gcore::System *d_sys;
    int d_mol;
    int d_atom;

    struct Exception: public gromos::Exception{
      Exception(const std::string &what):
	gromos::Exception("AtomSpecifier", what){}
    };

  };

  /**
   * Class SolventSpecAtom
   * purpose: interface to access a specific solvent atom (in an AtomSpecifier)
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class SolventSpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class SolventSpecAtom : public SpecAtom
  {
  public:
    SolventSpecAtom(gcore::System &sys, int m, int a) : SpecAtom(sys, m, a) {}
    SolventSpecAtom(SolventSpecAtom const &s) : SpecAtom(s) {}

    virtual ~SolventSpecAtom() {};

    virtual spec_type type() const { return spec_solvent; }

    virtual SpecAtom * clone() { return new SolventSpecAtom(*this); }

    virtual int mol() const { return d_mol; }
    virtual int atom()const { return d_atom; }

    virtual gmath::Vec & pos() {
      if(d_atom > d_sys->sol(0).numPos())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).pos(d_atom);
    }
    virtual gmath::Vec const & pos() const {
      if(d_atom > d_sys->sol(0).numPos())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).pos(d_atom);
    }
    virtual gmath::Vec & cosDisplacement() {
      if(d_atom > d_sys->sol(0).numPos())
  throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).cosDisplacement(d_atom);
    }
    virtual gmath::Vec const & cosDisplacement() const {
      if(d_atom > d_sys->sol(0).numPos())
  throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).cosDisplacement(d_atom);
    }
    virtual gmath::Vec & vel() {
      if(d_atom > d_sys->sol(0).numVel())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).vel(d_atom);
    }
    virtual gmath::Vec const & vel() const {
      if(d_atom > d_sys->sol(0).numVel())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).vel(d_atom);
    }

    virtual std::string name()const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).name();
    }

    virtual int iac() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).iac();
    }
    virtual double charge() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).charge();
    }
    virtual bool isPolarisable() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).isPolarisable();
    }
    virtual double poloffsiteGamma() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).poloffsiteGamma();
    }
    virtual int poloffsiteI() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).poloffsiteI();
    }
    virtual int poloffsiteJ() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).poloffsiteJ();
    }
    virtual double cosCharge() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).cosCharge();
    }
    virtual double mass() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).mass();
    }
    virtual int resnum() const
    {
      return d_atom / d_sys->sol(0).topology().numAtoms();
    }
    virtual std::string resname() const
    {
      return "SLV";
    }

  private:

  };

  /**
   * Class VirtualSpecAtom
   * purpose: interface to access a specific atom (in an AtomSpecifier)
   * unified access to solute, solvent and virtual atoms.
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class VirtualSpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class VirtualSpecAtom : public SpecAtom
  {
  public:
    VirtualSpecAtom(gcore::System &sys, std::string s, VirtualAtom::virtual_type t);
    VirtualSpecAtom(gcore::System &sys, std::string s, int x, VirtualAtom::virtual_type t);
    VirtualSpecAtom(VirtualSpecAtom const &s) : SpecAtom(s), d_va(s.d_va), d_pos(s.d_pos) {}

    virtual ~VirtualSpecAtom() {};

    virtual spec_type type() const { return spec_virtual; }

    virtual SpecAtom * clone() { return new VirtualSpecAtom(*this); }

    virtual int mol() const { /* throw Exception(" accessing VA mol");*/ return d_mol; }
    virtual int atom()const { /* throw Exception(" accessing VA atom"); */ return d_atom; }

    virtual gmath::Vec & pos() {return d_pos = d_va.pos(); }
    virtual gmath::Vec const & pos() const { return d_pos = d_va.pos(); }
    virtual gmath::Vec & vel() { throw Exception(" accessing VA velocitiy");}
    virtual gmath::Vec const & vel() const { throw Exception(" accessing VA velocitiy"); }

    virtual std::string name()const
    {
      return "VA";
    }

    virtual int iac() const
    {
      throw Exception(" accessing VA iac");
      // return 0;
    }
    virtual double charge() const
    {
      throw Exception(" accessing VA charge");
      // return 0;
    }
    virtual double mass() const
    {
      throw Exception(" accessing VA mass");
      // return 0.0;
    }
    virtual int resnum() const
    {
      throw Exception(" accessing VA resnum");
      // return 0;
    }
    virtual std::string resname() const
    {
      throw Exception(" accessing VA resname");
      // return "VA";
    }

    virtual void setSystem(gcore::System &sys);

    virtual std::string toString() const
    {
      return d_va.toString();
    }
    virtual utils::AtomSpecifier conf();

    virtual int virtualType()
    {
      return d_va.type();
    }

  protected:
    VirtualAtom d_va;
    mutable gmath::Vec d_pos;
  };

  /**
   * @class AtomSpecifier
   * @author C. Oostenbrink, M. Christen
   * @ingroup utils
   *
   * Class AtomSpecifier
   * purpose: contains specific atoms of the system, keeping track of
   * molecule and atom numbers
   *
   * Description:
   * This class defines (and implements) a general form to
   * access atoms in a system.
   *
   * @section AtomSpecifier Atom Specifier
   * The AtomSpecifier can be used to look over a specific set of atoms,
   * possibly spanning different molecules. Please be aware that some shells
   * may modify the brackets and other special characters. Therefore, quotes
   * should be used when using atom specifiers as command line arguments.
   * A 'specifier' is a string with following format:
   *
   * @subsection atomspec_mol_atom Molecules and Atoms
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <mol>[-<mol>]:<atom>[-<atom>] @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<mol\></span> is a
   *   molecule number, or 'a' for all molecules
   * - <span style="color:darkred;font-family:monospace">\<atom\></span> is an
   *   atom number, or 'a' for all atoms, or an atom name
   *
   * Molecules and atoms are seperated by a colon (:). Ranges of molecules, atoms
   * are seperated by a dash (-). Sets of atoms are seperated by a comma (,).
   * Molecules are specified by their molecule number or "a" or "s" which is a
   * shortcut notation for all molecules or solvent molecules respectively.
   * Atoms are specified by either their names, atom numbers or "a" which is a
   * shortcut notation for all all atoms.
   *
   * For example:
   * - @verbatim 2:5 @endverbatim means atom 5 of molecule 2.
   * - @verbatim 1:3-5 @endverbatim means atoms 3, 4 and 5 of molecule 1.
   * - @verbatim 2:3,5,8 @endverbatim means atoms 3, 5 and 8 of molecule 2.
   * - @verbatim 3:a @endverbatim means all atoms of molecule 3.
   * - @verbatim a:CA @endverbatim means all CA atoms of all molecules.
   * - @verbatim 1:CA,N,C @endverbatim means all CA, N or C atoms of molecule 1.
   *
   * @subsection atomspec_res Residues
   *
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <mol>[-<mol>]:res(<residue>:<atom>[,<atom>...]) @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<residue\></span> is a
   *   residue number or a residue name
   *
   * By default the residue information is ignored and the atom specifier only
   * consists of molecule and atom information. However, the residue number or
   * name can by included in the atom specifier by using the res() directive.
   * Within the res() directive the residue and atom information are sperated by
   * a colon (:). The residue information may consist of sets of residue names
   * or ranges and sets of residue numbers. Please note that only sets of atoms
   * (seperated by a comma) and no ranges of atoms can be used when using
   * res().
   *
   * For example:
   * - @verbatim 1:res(5:C) @endverbatim means atom C from residue 5 of molecule 1.
   * - @verbatim 1:res(1-5:CA) @endverbatim means all CA atoms from residues 1 to 5 of molecule 1.
   * - @verbatim 1:res(3,5:1) @endverbatim means atom 1 of residue 3 or 5 of molecule 1.
   * - @verbatim 1:res(ILE:a) @endverbatim means all atoms of all ILE residues of molecule 1.
   * - @verbatim 1:res(SER,THR:C,N,CA) @endverbatim means all C, N or CA atoms of residues named
   * SER or THR of molecule 1.
   *
   * @subsection atomspec_va Virtual Atoms
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim va(<type>, <atomspec>) @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<type\></span> is a
   *   @ref VirtualAtom type.
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is
   *   an @ref AtomSpecifier
   *
   * Atom specifiers can make use of virtual atoms. The type of the virtual atom
   * is a type number or keyword and documented @ref VirtualAtom "here".
   * Depending on the virtual atom type the
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span> atom
   * specifier has be of a certain size.
   *
   * For example:
   * - @verbatim va(com,1:a) @endverbatim means centre of mass of molecule 1.
   * - @verbatim va(1,1:res(1:CA,N,CB,C)) @endverbatim means the alpha hydrogen
   *  of residue 1 in a protein.
   *
   * @subsection atomspec_file Reading From Files
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim file(<filename>) @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<filename\></span> is
   *   the output of the atominfo program
   *
   * Atom specifiers can be read from a file. Complex atom specifiers can be
   * generated by the program atominfo. Its output is saved and can be
   * read by using the <span style="color:darkred;font-family:monospace">file()
   * </span> atom specifier.
   *
   * For example:
   *  - @verbatim $ atominfo @topo ex.topo @atomspec 1:CA > ca.spec @endverbatim
   * @verbatim file(ca.spec)@endverbatim means all CA atoms of the first molecule
   *
   * @subsection atomspec_multi Multiple Atom Specifiers
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <atomspec>[;<atomspec>] @endverbatim
   * </b></span>
   * <br>
   * where
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   *
   * Multiple atom specifiers are seperated by a semicolon (;)
   *
   * For example:
   * - @verbatim 1:1;2:1 @endverbatim means the first atom of the first molecule
   * and the first atom of the second molecule.
   *
   * @subsection atomspec_excl Exclusion Atom Specifiers
   * - @verbatim not(<atomspec>) @endverbatim
   * - @verbatim minus(<atomspec>) @endverbatim
   *
   * <br>
   * where
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> is an
   *   @ref AtomSpecifier
   *
   * The atoms given using the keyword <span style="color:darkred;font-family:monospace">not</span>
   * are never included int the resulting atom specifier even if they are added
   * later. The atoms given using the keyword <span style="color:darkred;font-family:monospace">minus</span>
   * are removed form the current atom specifier but they can be added again in additional atom specifiers.
   * <span style="color:darkred;font-family:monospace">minus</span> cannot be used for solvent.
   *
   * For example:
   * - @verbatim 1:a minus(1:res(2:a)) @endverbatim means all atoms of the first
   *   molecule but without the second residue.
   * - @verbatim not(a:C) 1:a @endverbatim means all atoms of the first molecule
   *   but without any atom named "C".
   * - @verbatim 1:a minus(1:res(2:a)) 1:res(2:C) @endverbatim  means all atoms
   *   of the first molecule but without the second residue with the exception of the "C" atom.
   *
   * @subsection atomspec_no Empty Atom Specifiers
   * - @verbatim no @endverbatim
   *
   * In some cases an AtomSpecifier is required by a program but one does not want to include
   * any atoms in the computations. In this case, the keyword <span style="color:darkred;font-family:monospace">no</span>
   * is given. It stands for an empty set of atoms.
   *
   * For example:
   * - @verbatim no @endverbatim means no atoms at all
   * - @verbatim 1:1 no 2:1 @endverbatim means the first atom of the first and second molecule.
   *
   * <b>See also</b> @ref PropertySpecifier "Property Specifier"
   * @ref VectorSpecifier "Vector Specifier"
   *
   */
  class AtomSpecifier;
  class AtomSpecifier{
    mutable std::vector<int> d_solventType; //chris said mutable was ok ;)

    mutable std::vector<SpecAtom *> d_specatom;
    mutable AtomSpecifier * d_not_atoms;

    gcore::System *d_sys;
    mutable int d_nsm;

  public:
    // Constructors
    /**
     * AtomSpecifier standard constructor
     */
    AtomSpecifier();
    /**
     * AtomSpecifier Constructor
     * @param sys The AtomSpecifier needs to know about the system. It
     *            does not know about any atoms yet.
     */
    AtomSpecifier(gcore::System &sys);
    /**
     * AtomSpecifier Constructor
     * @param sys The AtomSpecifier needs to know about the system.
     * @param s   A string of the correct format. Usually this is provided
     *            by the user, so it is assumed to start numbering at 1
     */
    AtomSpecifier(gcore::System &sys, std::string s);
    /**
     * AtomSpecifier Constructor
     * @param sys The AtomSpecifier needs to know about the system.
     * @param s   A string of the correct format. Usually this is provided
     *            by the user, so it is assumed to start numbering at 1
     * @param x   substitution value
     */
    AtomSpecifier(gcore::System &sys, std::string s, int x);

    /**
     * copy constructor!
     */
    AtomSpecifier(AtomSpecifier const & as);
    /**
     * AtomSpecifier Destructor
     */
    ~AtomSpecifier();

    /**
     * Method to set the system the atom specifier is referring to
     * @param sys the system
     */
    void setSystem(gcore::System &sys);
    /**
     * Method to add parse a string to the AtomSpecifier.
     * @param s Is assumed to be user-specified, with numbering starting at 1
     */
    int addSpecifier(std::string s, int x=-1);
    /**
     * Method to add parse a string to the AtomSpecifier
     * without redundancy checks.
     * @param s Is assumed to be user-specified, with numbering starting at 1
     */
    int addSpecifierStrict(std::string s, int x=-1);
    /**
     * Method to add parse a string to the AtomSpecifier
     * without redundancy checks.
     * @param s Is assumed to be user-specified, with numbering starting at 1
     */
    int appendSpecifier(std::string s, int x=-1);
    /**
     * Method to add a single molecule to the AtomSpecifier
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule.
     */
    int addAtom(int m, int a);
    /**
     * Method to add a single molecule to the AtomSpecifier
     * without redundancy checks
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule.
     */
    int addAtomStrict(int m, int a);
    /**
     * Method to add an atom from a Gromos96 number.
     * Numbering is here assumed to be Gromos96 (starting at 1)
     * @param a Gromos96 atom number
     */
    int addGromosAtom(int a);
    /**
     * Method to add a complete molecule.
     * Numbering is here assumed to be gromos++ (starting at 0)
     * @param m Molecule number.
     */
    int addMolecule(int m);
    /**
     * Method to remove an atom from the AtomSpecifier.
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule
     */
    int removeAtom(int m, int a);
    /**
     * Method to remove an atom from the AtomSpecifier.
     *
     * @param i remove the atoms with index 1 in the specifier
     */
    int removeAtom(int i);
    /**
     * Method to find the index of a specific atom in the AtomSpecifier
     *
     * Numbering is assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule
     */
    int findAtom(int m, int a)const;
    /**
     * Method to add all atoms of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     */
    int addType(int m, std::string s);
    /**
     * Method to add all atoms of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     */
    int addTypeStrict(int m, std::string s);
    /**
     * Method to add atoms of the specified range
     * of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     * @param beg begin of range
     * @param end end of range
     */
    int addType(int m, std::string s, int beg, int end);
    /**
     * Method to add atoms of the specified range
     * of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     * @param beg begin of range
     * @param end end of range
     */
    int addTypeStrict(int m, std::string s, int beg, int end);
    /**
     * Method to add all atoms (in all molecules) with a certain name to the
     * AtomSpecifier
     * @param s Atom name that is to be added (e.g. CA)
     */
    int addType(std::string s);
    /**
     * Method to add solvent atoms of a type
     * @param s Atom name that is to be added
     */
    int addSolventType(std::string s);

    /**
     * Method to add solvent, as it was specified in the topology
     */
    int addSolventType() const;
    /**
     * Method to sort the atoms ascending order. Some applications might
     * need the atoms to be ordered. This is just a simple bubble sort
     */
    void sort();
    /**
     * Member operator = copies one AtomSpecifier into the other
     */
    AtomSpecifier &operator=(const AtomSpecifier &as);
    /**
     * Member operator + adds two AtomSpecifiers. Atoms that appeared in
     * both are only listed once.
     */
    AtomSpecifier operator+(const AtomSpecifier &as);

    /**
     * Accessor, returns the molecule number of the i-th atom in the
     * AtomSpecifier
     */
    int mol(int i)const;
    /**
     * Accessor, returns the atom number of the i-th atom in the
     * AtomSpecifier
     */
    int atom(int i)const;
    /**
     * Accessor, returns the gromos atom number of the i-th atom in the
     * AtomSpecifier
     */
    int gromosAtom(int i)const;
    /**
     * Accessor, returns a pointer to the vector containing the atom
     * numbers.
     */
    std::vector<SpecAtom *> & atom();
    /**
     * const Accessor, returns a pointer to the vector containing the atom
     * numbers.
     */
    std::vector<SpecAtom *> const & atom()const;
    /**
     * Accessor, returns the number of atoms in the AtomSpecifier
     */
    unsigned int size()const;
    /**
     * Accessor, returns if empty
     */
    bool empty()const;
    /**
     * Function to empty the AtomSpecifier
     */
    void clear();
    /**
     * Accessor, returns a pointer to the system on which the AtomSpecifier
     * is based
     */
    gcore::System *sys()const;
    /**
     * Accessor, returns a pointer to the coordinates of the i-th
     * atom in the AtomSpecifier
     */
    gmath::Vec *coord(int i);
    /**
     * Accessor, returns the coordinates of the i-th
     * atom in the AtomSpecifier
     */
    gmath::Vec & pos(int i);
    /**
     * Accessor, returns the coordinates of the i-th
     * atom in the AtomSpecifier
     * const accessor
     */
    gmath::Vec const & pos(int i)const;
    /**
     * Accessor, returns the coordinates of the cos displacement vector
     * for the i-th atom in the AtomSpecifier
     */
    gmath::Vec & cosDisplacement(int i);
    /**
     * Accessor, returns the coordinates of the cos displacement vector
     * for the i-th atom in the AtomSpecifier
     * const accessor
     */
    gmath::Vec const & cosDisplacement(int i)const;
    /**
     * Accessor, returns the velocities of the i-th
     * atom in the AtomSpecifier
     */
    gmath::Vec & vel(int i);
    /**
     * Accessor, returns the velocities of the i-th
     * atom in the AtomSpecifier
     * const accessor
     */
    gmath::Vec const & vel(int i)const;
    /**
     * Accesor, returns the atom name of the i-th atom in the AtomSpecifier
     */
    std::string name(int i)const;
    /**
     * Accessor, returns the Iac of the i-th atom in the AtomSpecifier
     */
    int iac(int i)const;
    /**
     * Accessotr, returns the vdW radius of the i-th atom in the AtomSpecifier
     */
    double radius(int i)const;
    /**
     * Accessor, returns the charge of the i-th atom in the AtomSpecifier
     */
    double charge(int i)const;
    /**
     * Accessor, returns wether the i-th atom in the AtomSpecifier is polarisable
     */
    bool isPolarisable(int i)const;
    /**
     * Accessor, returns the polarisation gamma of the i-th atom in the AtomSpecifier
     */
    double poloffsiteGamma(int i)const;
    /**
     * Accessor, returns the first atom for the off-site polarisation of the i-th atom in the AtomSpecifier
     */
    int poloffsiteI(int i)const;
    /**
     * Accessor, returns the second atom for the off-site polarisation of the i-th atom in the AtomSpecifier
     */
    int poloffsiteJ(int i)const;
    /**
     * Accessor, returns the cos charge of the i-th atom in the AtomSpecifier
     */
    double cosCharge(int i)const;
    /**
     * Accessor, returns the mass of the i-th atom in the AtomSpecifier
     */
    double mass(int i)const;
    /**
     * Accessor, returns the residue number of atom i in the AtomSpecifier
     */
    int resnum(int i)const;
    /**
     * Accessor, returns the residue name of the i-th atom
     * in the AtomSpecifier
     */
    std::string resname(int i)const;
    /**
     * Method, returns a vector of strings that would reproduce the
     * AtomSpecifier
     */
    std::vector<std::string> toString()const;
    /**
     * Method, returns a string of just the i-th atom in the
     * AtomSpecifier
     */
    std::string toString(int i)const;
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says AtomSpecifier, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what):
	gromos::Exception("AtomSpecifier", what){}
    };
    /**
     * Method, returns the number of atoms a solvent molecule has. This works before any solvent is read in from a file.
     */
    int numSolventAtoms() const;

  protected:
    //Internal function
    /**
     * Parse the arguments string into the AtomSpecifier
     */
    void parse(std::string s, int x=-1);
    /** Parse the arguments string into the AtomSpecifier without redundancy checks
     */
    void parseStrict(std::string s, int x=-1);
    /**
     * parse a single specifier (no ';')
     */
    void parse_single(std::string s, int x=-1);
    /**
     * parse a single specifier (no ';')
     */
    void parse_singleStrict(std::string s, int x=-1);
    /**
     * parse the molecules, deliver in mol
     */
    void parse_molecule(std::string s, std::vector<int> & mol, int x=-1);
    /**
     * parse the molecules, deliver in mol
     */
    void parse_moleculeStrict(std::string s, std::vector<int> & mol, int x=-1);
    /**
     * parse the atom part
     * can also include residues
     */
    void parse_atom(int mol, std::string s, int x=-1);
    /**
     * parse the atom part
     * can also include residues
     */
    void parse_atomStrict(int mol, std::string s, int x=-1);
    /**
     * parse an atom from a range of possible atoms
     * the range is either all atoms of the molecule mol
     * or just the atoms inside a residue of the molecule mol
     */
    void parse_atom_range(int mol, int beg, int end, std::string s, int x=-1);
    /**
     * parse an atom from a range of possible atoms
     * the range is either all atoms of the molecule mol
     * or just the atoms inside a residue of the molecule mol
     */
    void parse_atom_rangeStrict(int mol, int beg, int end, std::string s, int x=-1);
    /**
     * parse virtual atoms.
     * standard virtual types plus
     * centre of geometry (cog) and centre of mass (com)
     */
    void parse_va(std::string s, int x = -1);
    /**
     * parse atom removal (minus(spec))
     */
    void parse_minus(std::string s, int x = -1);
    /**
     * parse atominfo file.
     */
    void parse_atominfo(std::string s);
    /**
     * parse atominfo file without redundancy checks.
     */
    void parse_atominfoStrict(std::string s);
    /**
     * parse residues
     */
    void parse_res(int mol, std::string res, std::string atom, int x=-1);
    /**
     * parse residues
     */
    void parse_resStrict(int mol, std::string res, std::string atom, int x=-1);
    /**
     * parse residues if specified by type
     */
    void parse_res_type(int mol, std::string res, std::string s, int x=-1);
    /**
     * Adds an atom to the atom specifier, checks whether it is in already
     */
    void _appendAtom(int m, int a);
    /**
     * Adds an atom to the atom specifier, checks whether it is in already
     * only used in solvent expansion, therefore const
     */
    void _appendSolvent(int m, int a)const;
    /**
     * Compares the string s to the name of the a-th atom in the m-th
     * molecule of the system. If they are the same, it returns true. Only
     * the characters before a possible '?' in s are compared.
     */
    bool _checkName(int m, int a, std::string s)const;
    /**
     * Compares the string s to the name of the residue
     * the  a-th atom in the m-th molecule of the system belong to.
     * If they are the same, it returns true. Only
     * the characters before a possible '?' in s are compared.
     */
    bool _checkResName(int m, int a, std::string s)const;
    /**
     * Method that expands specified solvent types to the real number of
     * solvent molecules. Is called from any accessor number of solvent
     * molecules has changed
     */
    int _expandSolvent()const;
    /**
     * Method to remove an atom from the AtomSpecifier.
     * only used to remove solvents before expanding
     * therefore const!
     * @param i remove the atoms with index 1 in the specifier
     */
    int _unexpandSolvent(int i)const;
    /**
     * Tells you if the number of solvent molecules in the system has changed
     * if so, the accessors should re-expand the Solvent Types.
     */
    bool _expand()const;
    /**
     * special function for sorting, to replace 'larger than' in comparisons
     * makes sure that solvent atoms will be last
     * @param i index of that is compared to
     * @param m molecule number
     * @param a atom number
     * @return bool returns true if atom i should come after m:a
     *                      false if atom i should come before m:a
     */
    bool _compare(int i, int m, int a)const;
    /**
     * find the atoms belonging to residue res of molecule mol
     */
    void res_range(int mol, int res, int &beg, int &end);

};
  //inline functions and methods

  inline std::vector<SpecAtom *> & AtomSpecifier::atom()
  {
    return d_specatom;
  }

  inline std::vector<SpecAtom *> const & AtomSpecifier::atom()const
  {
    return d_specatom;
  }

  inline gcore::System *AtomSpecifier::sys()const
  {
    return d_sys;
  }

  inline VirtualSpecAtom::VirtualSpecAtom(gcore::System &sys,
					  std::string s,
					  VirtualAtom::virtual_type t)
    : SpecAtom(sys, -3, -1),
      d_va(sys, AtomSpecifier(sys, s), t)
  {
  }

  inline VirtualSpecAtom::VirtualSpecAtom(gcore::System &sys,
					  std::string s, int x,
					  VirtualAtom::virtual_type t)
    : SpecAtom(sys, -3, -1),
      d_va(sys, AtomSpecifier(sys, s, x), t)
  {
  }

  inline void VirtualSpecAtom::setSystem(gcore::System &sys)
  {
    d_va.setSystem(sys);
    SpecAtom::setSystem(sys);
  }

  inline AtomSpecifier VirtualSpecAtom::conf()
  {
    return d_va.conf();
  }
  inline AtomSpecifier SpecAtom::conf()
  {
    AtomSpecifier as;
    return as;
  }

  inline int AtomSpecifier::numSolventAtoms() const{
    return d_solventType.size(); //this also works before any coordinates are read in
  }
}

#endif

