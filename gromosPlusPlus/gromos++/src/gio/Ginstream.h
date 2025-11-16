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
 * @file Ginstream.h
 * basic input stream class.
 */

#ifndef INCLUDED_GINSTREAM_H
#define INCLUDED_GINSTREAM_H

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

namespace gio {

  /**
   * @class Ginstream
   * @ingroup gio
   * the basic block input stream
   */
  class Ginstream {

  public:
    /**
     * Default Constructor
     */
    Ginstream():_is(){ _has_version=false;};
    
    
    /**
     * Constructor with an existing stream
     */
    Ginstream(std::ifstream& is) { stream(is); _name=""; _has_version=false;} 
    /**
     * Copy constructor
     */
    Ginstream(Ginstream &gin);
    
    /**
     * Constructor that tries to open a file
     */
    Ginstream(const std::string &s, std::ios::openmode mode = std::ios::in);
    

    /*
     * Accessors to the input stream
     */
    std::istream& stream() { return *_is; }
    void stream(std::istream& is, bool ver = false) { 
        _is = &is; 
        readTitle(); 
        if (ver) {
            readVersion();
        }
    }
    void open(const std::string s, std::ios::openmode mode = std::ios::in);
    void close();
    
    /**
     * Read a title block from the input stream,
     * concatenate and store it in title.
     */
    void readTitle();

    /**
     * Accessor, returns the title of the Ginstream
     */
    std::string title();
    
    /**
     * Read a version block from the input stream,
     * concatenate and store it in version.
     */
    void readVersion();

    /**
     * Returns true if the Ginstream has a version
     */
    bool has_version() const ;

    /**
     * Accessor, returns the version of the Ginstream
     */
    std::string version() const ;

    /**
     * Accessor, returns the name of the file (if constructed or
     * opened with a filename)
     */
    std::string name() const ;
    
    /** 
     * The function gio::Ginstream::getline provides an override of 
     * std::getline in that it retrieves the next line (separated 
     * by const char& sep) from the input stream, which, after 
     * having been stripped of comments (indicated by 
     * const char& comm), is not empty.
     */
    std::istream& getline(std::string& s, 
			  const char& sep = '\n',
			  const char& comm = '#');
    
    /**
     * The function gio::Ginstream::getblock retrieves the next block  
     * (separated by const string& sep) using io::getline.
     *
     * It throws a runtime_error if the stream's good bit is unset 
     * before it's finished reading a block.
     *
     * Finally, the vector<string> it writes to is resized to the 
     * number of strings read.
     */
    std::istream& getblock(std::vector<std::string>& b, 
			   const std::string& sep = "END");

    /**
     * The function gio::Ginstream::skipblock retrieves the next block  
     * (separated by const string& sep) using io::getline.
     *
     * It throws a runtime_error if the stream's good bit is unset 
     * before it's finished reading a block.
     *
     * the block is discarded while read.
     */
    std::istream& skipblock(const std::string& sep = "END");

  protected:
    std::istringstream _lineStream;

  private:
    std::istream* _is;
    std::string _title;
    std::string _name;
    bool _has_version;
    std::string _version;
  };    

  /**
   * The gio::concatenate utility allows the concatenation
   * of entries in a vector<string> into just one string. The
   * resulting entries are separated by const char& sep 
   * (defaulting to a newline character).
   */
  std::string& concatenate(std::vector<std::string>::const_iterator begin,
			   std::vector<std::string>::const_iterator end,
			   std::string& s,
			   const char& sep = '\n');
} // gio

#endif
