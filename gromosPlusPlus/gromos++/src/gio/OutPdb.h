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

// gio_OutPdb.h

#ifndef INCLUDED_GIO_OUTPdb
#define INCLUDED_GIO_OUTPdb

#include <string>

#include "OutCoordinates.h"

namespace gcore{
  class System;
}

namespace gio{
  class OutPdb_i;
  /**
   * Class OutPdb
   * is of type OutCoordinates and defines how a PDB-file is written out
   * 
   * @class OutPdb
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @ingroup gio
   */
  class OutPdb: public OutCoordinates{
    OutPdb_i *d_this;
    // prevent copying and assignment
    OutPdb(const OutPdb &);
    OutPdb &operator=(const OutPdb&);
    std::string flavour;
    double factor;
    bool renumber;
  public:
    OutPdb(std::string flavour = "pdb", double factor = 10.0, bool renumber=false);
    OutPdb(std::ostream &os, std::string flavour = "pdb", double factor = 10.0, bool renumber=false);
    ~OutPdb();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    OutPdb &operator<<(const gcore::System &sys);
    OutPdb &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
