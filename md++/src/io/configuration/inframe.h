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
 * @file inframe.h
 * frame based input / output for trajectories
 */

#ifndef INCLUDED_INFRAME_H
#define INCLUDED_INFRAME_H

#include "../instream.h"

namespace io {

  class In_Frame : public GInStream {

  public:

    /*
     * default constructor
     */
    In_Frame() : GInStream(){} 

    /*
     * Constructor
     */
    In_Frame(std::istream& is) : GInStream(is) {
      io::getline(*_is, _next_block);
      _first_block = _next_block;
    } 
    
    /*
     * Accessors to the input stream
     */
    void stream(std::istream& is) { 
      GInStream::stream(is);io::getline(*_is, _next_block);
      _first_block = _next_block;  
    }

    /**
     * read frame.
     */
    int read_frame();
    

  protected:
    /**
     * next block name in file.
     */
    std::string _next_block;
    /**
     * first block name in file (after title block)
     */
    std::string _first_block;
    
  };

} // io

#endif
