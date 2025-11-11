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

// gio_OutG96.h

#ifndef INCLUDED_GIO_OUTG96
#define INCLUDED_GIO_OUTG96

#include <string>

#include "OutCoordinates.h"

namespace gcore{
  class System;
  class Box;
}

namespace utils {
  class AtomSpecifier;
}

namespace gio{
  class OutG96_i;
  /**
   * Class OutG96
   * is of type OutCoordinates and defines how a gromos96 trajectory is 
   * printed out (POSITIONRED block) 
   *
   * @class OutG96
   * @author R. Buergi
   * @ingroup gio
   */
  class OutG96: public OutCoordinates{
    OutG96_i *d_this;
    // prevent copying and assignment
    OutG96(const OutG96 &);
    OutG96 &operator=(const OutG96&);
  public:
    OutG96();
    OutG96(std::ostream &os);
    ~OutG96();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    void writeGenBox(const gcore::Box &box);
    void writeTriclinicBox(const gcore::Box &box);
    OutG96 &operator<<(const gcore::System &sys);
    OutG96 &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
