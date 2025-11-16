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
 * @file Ginstream.cc
 * basic input stream class definition.
 */

#include "Ginstream.h"

#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>    // std::remove_if
#include <cstdio>

#include "gzstream.h"
#include "../gromos/Exception.h"

template<class size_type>
inline std::basic_string<size_type>&
trim_right( std::basic_string<size_type>& str )
{
    return( str = str.substr( 0, str.find_last_not_of( ' ' ) + 1 ) );
}


template<class size_type>
inline std::basic_string<size_type>&
trim( std::basic_string<size_type>& str )
{
  if (str.find_first_not_of( ' ' ) == std::string::npos) return (str = "");  // line of whitespaces
  return( trim_right( str ) );  //remove trailing whitespace from string
}


gio::Ginstream::Ginstream(const std::string &s, std::ios::openmode mode)
  :_is(NULL)
{
  open(s, mode);
  _has_version = false;
}

gio::Ginstream::Ginstream(gio::Ginstream &gin)
{
  _is = gin._is;
  _title = gin._title;
  _name = gin._name;
  _has_version = gin._has_version;
  _version = gin._version;
}

void gio::Ginstream::open(const std::string s, std::ios::openmode mode)
{
  igzstream *gis = new igzstream(s.c_str(), mode);
  if (gis == NULL){
    throw gromos::Exception("Ginstream", "could not create a std::ifstream( " + s + ")");
  }
  if(!gis->good()){
    throw gromos::Exception("Ginstream", "Could not open file '" + s + "'");
  }  
  if(!gis->is_open()){
    throw gromos::Exception("Ginstream", "could not open file '" + s + "'");
  }
  
  stream(*gis, true);
  if (!_has_version) {
    // We tried to read ENEVERSION, but apparently, there was none.
    // Our stream isn't usable anymore - we started reading the 
    // next block - probably TIMESTEP.
    // Rewind operations (seekg() and the like) don't seem to work - 
    // so let's recreate it. (Yeah, not THAT elegant...)
    gis->close();
    delete gis;
    igzstream *gis = new igzstream(s.c_str(), mode);
    if (gis == NULL){
      throw gromos::Exception("Ginstream", "could not create a std::ifstream( " + s + ")");
    }
    if(!gis->good()){
      throw gromos::Exception("Ginstream", "Could not open file '" + s + "'");
    }  
    if(!gis->is_open()){
      throw gromos::Exception("Ginstream", "could not open file '" + s + "'");
    }
    
    stream(*gis, false);
  }
  _name=s;
}


void gio::Ginstream::close()
{
  // std::cerr << "closing _is" << std::endl;

  std::ifstream * ifs = dynamic_cast<std::ifstream *>(_is);
  if (ifs != NULL)
    ifs->close();
  else{
    igzstream * igzs = dynamic_cast<igzstream *>(_is);
    igzs->close();
  }

  // _is->close();
  _title="";
  _name="";
  delete _is;
}

void gio::Ginstream::readTitle() {

  std::vector<std::string> _b;

  getblock(_b);
  std::string title = _b[0];
  title.resize (5);
  if (title != "TITLE")
    throw gromos::Exception("Ginstream",
                            "TITLE block expected. Found: " +  _b[0]);
  // Storing files on Windows causes weird end line characters
  if (_b[0] != "TITLE")
    throw gromos::Exception("Ginstream",
                            "ERROR Trailing data after keyword TITLE; possibly MS DOS newline character");
  _title = gio::concatenate(_b.begin() + 1, _b.end() - 1, _title);
}

std::string gio::Ginstream::name() const
{
  return _name;
}

std::string gio::Ginstream::title() 
{
  return _title;
}

std::istream& gio::Ginstream::getline(std::string& s, 
					     const char& sep,
					     const char& comm){
  //unsigned short int ii;
  std::string::size_type ii;
  
  while (_is->good()) {
    std::getline(*_is, s, sep); // the line is stored in s
    //ii = std::find(s.begin(), s.end(), comm) - s.begin();
    ii=s.find(comm,0); // ii = position where a comment character is found

    if(!s.size()) continue;                 // empty line. goto next line
    else if(ii == std::string::npos) break; // no comment. finished getting the line
    else if (!ii) continue;                 // comment on first position. goto next line
    else {
      s.erase(s.begin() + ii, s.end()); // delete comment from s
      if (!trim_right(s).size()) continue;  // line with comment only. there were only " " in the line. goto next line
      break;
    }
  }
  s = trim(s);  // remove trailing whitespace from line

  return *_is;
}

std::istream& gio::Ginstream::getblock(std::vector<std::string>& b, 
					      const std::string& sep)
{
  if (!b.size())
    b.push_back("");
  std::vector<std::string>::iterator dest = b.begin();

  while (1) {

    if (dest == b.end()) {
      b.push_back("");  //prepare the new vector item so it can receive a line
      dest = b.end() - 1;  // set the iterator to the last string
    }       
    getline(*dest); //read the line
   
    if(_is->eof()) {
      //solution for END<eof> (without newline after END):
      if(dest->std::string::empty())  // only rollback if the line is empty. getline will remove any whitespace if the line only contained whitespace. does not work for tab
        --dest;
      break;
    }
    
    if (dest->find(sep)==0){  // the line startswith separator END
      break;
    }
   
    if (!_is->good()) 
      throw gromos::Exception("getblock", "error reading block."+*b.begin());
        
    ++dest;
  }
  
  ++dest;
  b.erase(dest, b.end());

  return *_is;
}

std::istream& gio::Ginstream::skipblock(const std::string& sep)
{
  std::string s;

  while (1) {

    getline(s);

    if(_is->eof()) {
      break;
    }
    
    if (s.find(sep) == 0)
      break;

    if (!_is->good()) 
      throw gromos::Exception("getblock", "error skipping block." + s);
  }
  
  return *_is;
}


std::string& gio::concatenate(
		std::vector<std::string>::const_iterator begin,
		std::vector<std::string>::const_iterator end,
		std::string& s,
		const char& sep) {
  //s.clear();
  s="";
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }
  
  return s;
}

void gio::Ginstream::readVersion(){  
  std::vector<std::string> _b;

  getblock(_b);
  if (_b[0] != "ENEVERSION") {
    _has_version = false;
    return;
  }
  // else, save new version and set bool switch
  _version = gio::concatenate(_b.begin() + 1, _b.end() - 1, _version);
  // we're ignoring any whitespaces - less error-prone
  _version.erase( std::remove_if( _version.begin(), _version.end(), ::isspace ), _version.end() );

  _has_version = true;
}

bool gio::Ginstream::has_version() const {
  return _has_version;
}

std::string gio::Ginstream::version() const {
  if (_has_version) {
    return _version;
  } else {
    return "";
  }
}

