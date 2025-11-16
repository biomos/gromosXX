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
 * @file VectorSpecifier.h
 * VectorSpecifier methods
 */

// Class that contains a vector
// specified by cartesian coordinates,
// polar coordinates or
// 

#ifndef INCLUDED_UTILS_VECTORSPECIFIER
#define INCLUDED_UTILS_VECTORSPECIFIER

#include <vector>
#include <string>
#include <map>

// minimal complete headers
#include "../utils/AtomSpecifier.h"
#include "../gromos/Exception.h"

namespace gmath
{
  class Vec;
}

namespace bound
{
  class Boundary;
}

namespace utils
{
  class Value;
  
  /**
   * @class VectorSpecifier
   * @author C. Oostenbrink, M. Christen
   * @ingroup utils
   *
   * Class VectorSpecifier
   * purpose: specify vectors as absolute values or by giving
   * one or two atoms.
   *
   * @section VectorSpecifier Vector Specifier
   * There are three different ways of specifying a vector:
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim cart(<x>,<y>,<z>) @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<x\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<y\></span> and
   *   <span style="color:darkred;font-family:monospace">\<z\></span> are
   *   cartesian coordinates.
   *
   * This creates a vector @f$\vec{x}@f$ with the cartesian coordinates 
   * <span style="color:darkred;font-family:monospace">\<x\></span>, 
   * <span style="color:darkred;font-family:monospace">\<y\></span> and
   * <span style="color:darkred;font-family:monospace">\<z\></span>.
   * 
   * For example:
   * - @verbatim cart(2,5,1) @endverbatim means the vector @f$(2, 5, 1)@f$.
   *
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim polar(<r>,<alpha>,<beta>) @endverbatim
   * </b></span>
   * where
   * - <span style="color:darkred;font-family:monospace">\<r\></span> is the 
   *  length @f$r@f$ of the vector.
   * - <span style="color:darkred;font-family:monospace">\<alpha\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<beta\></span> are 
   *   angles @f$\alpha@f$ and @f$\beta@f$ in degrees.
   *
   * This creates a vector @f$\vec{x}@f$ with the cartesian coordinates
   * @f[\vec{x} = (r\cos(\alpha)\cos(\beta), r\sin(\alpha), -r\cos(\alpha)\sin(\beta)) @f]
   *
   * For example:
   * - @verbatim polar(2.5,45.0,90.0) @endverbatim means the vector @f$(0, 2.5, 0)@f$.
   *
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim atom(<atomspec>) @endverbatim
   * </b></span>
   * <br>
   * where
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   *
   * An atom specifier must contain one or two atoms, virtual atoms
   * are allowed. If only one atom is given, its position is taken as the 
   * vector. If two atoms are given, the vector pointing from atom one to
   * atom two is taken as the vector. The periodic boundary conditions are taken
   * into account (nearest image distance vector).
   *
   * For example:
   * - @verbatim atom(1:1) @endverbatim means the position vector of atom 1 of
   *   molecule 1.
   * - @verbatim atom(1:1,5) @endverbatim means the vector pointing from atom 1
   *   to atom 5 of molecule 1.
   * - @verbatim atom(va(cog,1:a);2:1) @endverbatim means the vector pointing 
   *   from the centre of geometry of molecule 1 to the first atom of molecule 2.
   * 
   * <b>See also</b> @ref AtomSpecifier "Atom Specifier"
   * @ref PropertySpecifier "Property Specifier"
   */
  class VectorSpecifier{

    AtomSpecifier d_atomspec;
    mutable gmath::Vec d_vec;
    bound::Boundary * d_pbc;
    
  public: 
    typedef std::map<std::string, utils::Value> var_type;

    // Constructors
    /** 
     * VectorSpecifier standard constructor
     */
    VectorSpecifier() : d_atomspec(), d_vec()
    {}

    /**
     * VectorSpecifier Constructor
     * @param sys The VectorSpecifier needs to know about the system. It 
     *            does not know about any atoms yet.
     */
    VectorSpecifier(gcore::System &sys, bound::Boundary * pbc)
      : d_atomspec(sys), d_vec(), d_pbc(pbc)
    {}
    
    /**
     * VectorSpecifier Constructor
     * @param sys The VectorSpecifier needs to know about the system.
     * @param pbc periodic boundary conditions
     * @param s   A string of the correct format. Usually this is provided
     *            by the user, so it is assumed to start numbering at 1
     * @param var variable substitutions while parsing
     */
    VectorSpecifier(gcore::System &sys, bound::Boundary * pbc, 
		    std::string s, var_type var = var_type());

    /**
     * copy constructor!
     */
    VectorSpecifier(VectorSpecifier const & vs);

    /**
     * VectorSpecifier Destructor
     */
    ~VectorSpecifier();

    /**
     * Method to set the system the atom specifier is referring to
     * @param sys the system
     */
    void setSystem(gcore::System &sys)
    {
      d_atomspec.setSystem(sys);
    }
    
    /**
     * set the boundary condition object.
     */
    void setBoundary(bound::Boundary *pbc)
    {
      d_pbc = pbc;
    }
    
    /**
     * Method to add parse a string to the VectorSpecifier.
     * @param s Is assumed to be user-specified, 
     * with numbering starting at 1
     */
    int setSpecifier(std::string s,
		     var_type var
		     = var_type());

    /**
     * Member operator = copies one VectorSpecifier into the other
     */
    VectorSpecifier &operator=(const VectorSpecifier &vs);

    /**
     * Accessor, returns the vector
     */    
    gmath::Vec const & operator()()const;
    /**
     * set vector to zero, empty atom specifier
     */
    void clear();
    /**
     * Accessor, returns a pointer to the system on which the VectorSpecifier
     * is based
     */
    gcore::System & sys();
    /**
     * boundary condition object accessor
     */
    bound::Boundary * pbc();
    
    /**
     * Method, returns a vector of strings that
     * would reproduce the
     * VectorSpecifier
     */
    std::string toString()const;
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says VectorSpecifier, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): 
	gromos::Exception("VectorSpecifier", what){}
    };
  protected:
    //Internal function
    /**
     * Parse the arguments string into the VectorSpecifier
     */
    void parse(std::string s, var_type & var);

    void parse_cartesian(std::string s);
    void parse_polar(std::string s);
    void parse_atom(std::string s, var_type & var);

  };
  
}

#endif

