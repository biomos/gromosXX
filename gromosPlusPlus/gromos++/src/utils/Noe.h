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

#ifndef INCLUDED_UTILS_NOE
#define INCLUDED_UTILS_NOE

#include <vector>

#include "../gromos/Exception.h"

namespace gmath {
  class Vec;
}

namespace gcore{
  class System;
}

namespace utils {
  class VirtualAtom;
}


namespace utils{
  class Noe_i;
  /**
   * Class Noe
   * The Noe class stores and analyses Noe information
   *
   * It stores the virtual atoms that define the NOE distance and the 
   * distances
   *
   * @class Noe
   * @author R. Buergi and M.A. Kastenholz
   * @ingroup utils
   */  
  class Noe{

    Noe_i *d_this;

    // not implemented
    Noe();
    Noe(const Noe &);
    Noe &operator=(const Noe&);
    
  
  public:
    /**
     * create an NOE from a line of GROMOS NOE calculation specification
     * @param sys the system to create the virtual atoms
     * @param line the line containing the distance restraint specification
     * @param dish carbon-hydrogen distance
     * @param disc carbon-carbon distance
     */
    Noe(gcore::System &sys, const std::string &line, double dish, double disc);
  
    /**
     * the distance corresponding to the NOE.
     * Periodic boundary conditions are not taken into account.
     */
    double distance(int i)const;
    /**
     * the distance vector corresponding to the NOE.
     * Periodic boundary conditions are not taken into account.
     */
    gmath::Vec distanceVec(int i)const;
    
    /**
     * get the virtual atoms
     * @param i the atom i occuring
     * @param ii of the distance ii
     */
    const utils::VirtualAtom & getAtom(int i, int ii) const;
   
    /**
     * the reference distance including the correction
     */
    double correctedReference(int i)const;

    /**
     * a string containing the GROMOS distance restraint specification line
     */
    std::string distRes(int i)const;

    /**
     * a string containing information about the NOE including
     * residue number, residue name, atom and molecule number and type.
     */
    std::string info(int i)const;


    /**
     * the number of distances
     */
    int numDistances()const;
    /**
     * the number of references
     */
    int numReferences()const;
    /**
     * the reference length of the NOE.
     * @param i index of the reference length (0)
     */
    double reference(int i)const;
    /**
     * the correction length for type
     * @param type the virtual atom type
     */
    double correction(int type);

    /**
     * set the correction length for a type
     * @param type virtual atom type
     * @param correction correction length
     */
    void setcorrection(int type, double correction);
    
    struct Exception: public gromos::Exception{
      Exception(const std::string &str): gromos::Exception("Noe", str){}
    };
    
  };

  /**
   * @class Noelib
   * A class to hold and parse NOE library information
   */
  class Noelib {
  public:
    std::string resname;
    std::string orgatomname;
    std::string gratomname;
    int NOETYPE;
    int NOESUBTYPE;

    Noelib(const std::string & A, const std::string & B, const std::string & C, 
    const std::string & D, const std::string & E = std::string("0")) {
      resname = A;
      orgatomname = B;
      gratomname = C;
      NOETYPE = atoi(D.c_str());
      NOESUBTYPE = atoi(E.c_str());
    }

    ~Noelib() {
    }
  };

  /**
   * parse an NOELIB block
   * @param buffer the buffer to parse
   * @param noelib the noe library
   */
  void parse_noelib(const std::vector<std::string> & buffer, std::vector<Noelib> & noelib);
  /**
   * create virtual atom(s) based on a type, subtype and atom number.
   * Automatically creates the r / l isomers.
   */
  std::vector<VirtualAtom*> getvirtual(int at, int type, int subtype, gcore::System &sys,
        double dish, double disc);

} /* namespace */

#endif
