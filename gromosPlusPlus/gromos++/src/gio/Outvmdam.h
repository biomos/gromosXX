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

// gio_Outvmdam.h

#ifndef INCLUDED_GIO_OUTvmdam
#define INCLUDED_GIO_OUTvmdam

#include <string>

#include "OutCoordinates.h"

namespace gcore{
  class System;
}

namespace gio{
  class Outvmdam_i;
  /**
   * Class Outvmdam
   * is of type OutCoordinates and defines how a trajectory should be
   * written out in "Amber Coordinates" that can be read by VMD.
   *
   * @class Outvmdam
   * @author M.K. Kastenholz
   * @ingroup gio
   */
  class Outvmdam: public OutCoordinates{
    Outvmdam_i *d_this;
    double factor;
    // prevent copying and assignment
    Outvmdam(const Outvmdam &);
    Outvmdam &operator=(const Outvmdam&);
  public:
    Outvmdam(double factor = 10.0);
    Outvmdam(std::ostream &os, double factor = 10.0);
    ~Outvmdam();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    Outvmdam &operator<<(const gcore::System & sys);
    Outvmdam &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
