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

// gio_OutGromacs.h

#ifndef INCLUDED_OUTGROMACS
#define INCLUDED_OUTGROMACS

#include<string>


namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  /**
   * Class OutGromacs
   * an outstream that writes a topology that looks similar to a 
   * gromacs topology
   * 
   * @class OutGromacs
   * @author B.C. Oostenbrink
   * @ingroup gio
   */
  class OutGromacs{
    std::string d_title;
    std::ostream &d_os;
    // prevent copying and assignment
    OutGromacs();
    OutGromacs(const OutGromacs&);
    OutGromacs &operator=(const OutGromacs&);
  public:
    /**
     * construct using an output stream
     */
    OutGromacs(std::ostream &os);
    ~OutGromacs();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the system and force-field parameters in gromacs format
     * @param sys the system
     * @param gff the force field
     */
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
  };
}

#endif
