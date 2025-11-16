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

// gio_OutTopology.h

#ifndef INCLUDED_OUTTOPOLOGY
#define INCLUDED_OUTTOPOLOGY

#include <string>
#include <set>

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  /**
   * Class OutTopology
   * an outstream that defines how a GROMOS topology should be written out
   *
   * @class OutTopology
   * @author B.C. Oostenbrink
   * @ingroup gio
   */
  class OutTopology{
    std::string d_title;
    std::ostream &d_os;
    
    // prevent copying and assignment
    OutTopology();
    OutTopology(const OutTopology&);
    OutTopology &operator=(const OutTopology&);
  public:
    /**
     * create from a stream
     */
    OutTopology(std::ostream &os);
    ~OutTopology();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the system and force field as a GROMOS topology to the stream
     * @param sys the system
     * @param gff the force field
     */
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
    /**
     * write the system and force field as a GROMOS96 topology to the stream
     * @param sys the system
     * @param gff the force field
     */
    void write96(const gcore::System &sys, const gcore::GromosForceField &gff);
  };
}
#endif
