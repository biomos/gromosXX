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
 * @file instream.h
 * basic input stream class.
 */

#ifndef INCLUDED_GINSTREAM_H
#define INCLUDED_GINSTREAM_H

namespace io {

  class GInStream {

  public:

    /*
     * default constructor
     */
    GInStream() : quiet(false), _auto_delete(false) {} 

    /*
     * Constructor
     */
    GInStream(std::istream& is) : quiet(false), _auto_delete(false) { stream(is); } 

    /**
     * Destructor.
     */
    ~GInStream() { if (_auto_delete) delete _is; }
    
    /*
     * Accessors to the input stream
     */
    std::istream& stream() { return *_is; }
    void stream(std::istream& is) { _is = &is; readTitle(); }

    /**
     * Read a title block from the input stream,
     * concatenate and store it in title.
     * @section title TITLE block
 @verbatim
TITLE
 your title
END
 @endverbatim
     */
    void readTitle();
    std::string title;
    /**
     * read the entire stream and store the blocks in the map.
     */
    void readStream();
    /**
     * auto delete accessor.
     */
    void auto_delete(bool b) { _auto_delete = b; }

    /**
     * output during reading?
     */
    bool quiet;

  protected:

    /**
     * temporary storage of individual block data while
     * the block is read
     */
    std::istringstream _lineStream;
    
    /**
     * stores the blocks if read_stream is called.
     */
    std::map<std::string, std::vector<std::string> > m_block;

    /**
     * the stream.
     */
    std::istream* _is;
    /**
     * delete the stream in the end?
     */
    bool _auto_delete;
    
  };

} // io

#endif
