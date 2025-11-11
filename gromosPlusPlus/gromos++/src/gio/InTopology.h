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

// gio_InTopology.h

#ifndef INCLUDED_GIO_INTOPOLOGY
#define INCLUDED_GIO_INTOPOLOGY

#include <string>

#include "../gromos/Exception.h"

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  class InTopology_i;
  /**
   * Class InTopology
   * defines an instream that can read in a GROMOS topology
   *
   * The data that is read in is split up into a System and a
   * GromosForceField
   *
   * @class InTopology
   * @ingroup gio
   * @author R. Buergi
   * @author B.C. Oostenbrink (massType, Solvent)
   * @sa gcore::System
   * @sa gcore::GromosForceField
   */
  class InTopology{
    InTopology_i *d_this;
    
  public:
    /**
     * open a topology file
     * @param file the topology file to open
     */
    InTopology(std::string file);
    ~InTopology();
    
    /**
     * access to the system that was read
     */
    const gcore::System &system()const;
    /**
     * access to the force field that was read
     */
    const gcore::GromosForceField &forceField()const;

    /**
     * access to the version string
     */
    const std::string &version()const;
    /**
     * access to the title
     */
    const std::string title()const;

    /**
     * The exception type for toplogy reading
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InTopology", what){}
    };
  };
}
#endif
