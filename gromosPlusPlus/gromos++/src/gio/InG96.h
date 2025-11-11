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

// gio_InG96.h

#ifndef INCLUDED_GIO_INCOORD
#define INCLUDED_GIO_INCOORD

#include <string>

#include "../gromos/Exception.h"

namespace gcore{
  class System;
}

namespace utils{
  class Time;
}

namespace gio{
  
  class InG96_i;
  /**
   * Class InG96
   * Defines an instream that can read any GROMOS coordinate file
   * or trajectory file.
   *
   * The instream can handle POSITION POSITIONRED VELOCITY and VELOCITYRED 
   * blocks in GROMOS files
   *
   * @class InG96
   * @ingroup gio
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @author M. Christen (velocities)
   * @sa gcore::System
   */
  class InG96{
    /**
     * pIMPL
     */
    InG96_i *d_this;
    /**
     * skip
     */
    int d_skip;
    /**
     * stride
     */
    int d_stride;
    /**
     * end of file during striding (or skipping)
     */
    int d_stride_eof;
    
    /**
     * copy constructor
     * not implemented
     */
    InG96(const InG96&);
    /**
     * operator =
     * not implemented
     */
    InG96 &operator=(const InG96&);
  public:
    /**
     * Constructor
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    InG96(int skip=0, int stride=1);
    /**
     * Constructor
     * @param name   : trajectory file
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    InG96(const std::string &name, int skip=0, int stride=1);
    /**
     * Destructor
     */
    ~InG96();

    // Methods
    /**
     * select SOLUTE, SOLVENT, ALL
     */
    void select(const std::string &thing);
    /**
     * open a trajectory file
     * skip and stride continue
     */
    void open(const std::string &name);
    /**
     * close
     */
    void close();
    /**
     * read a frame
     */
    InG96 &operator>>(gcore::System &sys);
    /** 
     * get the time information from a frame
     */
    InG96 &operator>>(utils::Time &time);

    // Accessors
    /**
     * title (of the trajectory)
     */
    std::string title()const;
    /**
     * name?
     */
    std::string name()const;
    /**
     * end of file
     */
    bool eof()const;
    
    /**
     * take every n-th frame
     * starts with frame 0, then n, 2n, 3n, ...
     */
    int stride()const;
    /**
     * skip n frames
     * after skipping, skip will be 0.
     * to skip again, just set it again...
     */
    int skip()const;
    
    /**
     * set stride
     */
    void stride(int stride);
    /**
     * set skip
     */
    void skip(int skip);

    /**
     * encountered end of file
     * during striding?
     */
    bool stride_eof()const;

    //Exceptions
    /**
     * Exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("InG96", what_arg){}
    };
  };
}
#endif
