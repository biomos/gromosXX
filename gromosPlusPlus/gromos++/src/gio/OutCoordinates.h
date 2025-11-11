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

// gio_OutCoordinates.h

#ifndef INCLUDED_GIO_OUTCOORDINATES
#define INCLUDED_GIO_OUTCOORDINATES

#include<string>

namespace gcore {
  class System;
}

namespace utils {
  class AtomSpecifier;
}

namespace gio{
  /**
   * Class OutCoordinates
   * defines some basic features for an output stream that writes out
   * GROMOS coordinate or trajectory files
   *
   * @class OutCoordinates
   * @author R. Buergi
   * @ingroup gio
   */
  class OutCoordinates{
    // prevent copying and assignment
    OutCoordinates(const OutCoordinates &);
    OutCoordinates &operator=(const OutCoordinates&);
  public:
    OutCoordinates(){}
    virtual ~OutCoordinates();
    /**
     * open an output stream
     */
    virtual void open(std::ostream &os)=0;
    /**
     * close the output stream
     */
    virtual void close()=0;
    /**
     * write the title string
     */
    virtual void writeTitle(const std::string &title)=0;
    /**
     * write the time and step information
     */
    virtual void writeTimestep(const int step, const double time)=0;
    /**
     * select only parts of the system
     * @param thing ALL for everything, SOLVENT for solvent only, or anything else for solute
     */
    virtual void select(const std::string &thing)=0;
    /**
     * write a system to the stream
     */
    virtual OutCoordinates &operator<<(const gcore::System & sys)=0;
    /**
     * write an AtomSpecifier to the stream
     */
    virtual OutCoordinates &operator<<(const utils::AtomSpecifier & atoms)=0;
  };
}

#endif
