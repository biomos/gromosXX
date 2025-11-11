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

// utils_RestrTraj.cc
#include "RestrTraj.h"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "groTime.h"
#include "../gio/Ginstream.h"

enum blocktype {titleblock, timestep, jvaluereseps, xrayrvalue};

typedef std::map<std::string,blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TITLE", titleblock),
                         BT("TIMESTEP", timestep),
                         BT("JVALUERESEPS", jvaluereseps),
                         BT("XRAYRVALUE", xrayrvalue)};


static std::map<std::string,blocktype> BLOCKTYPE(blocktypes,blocktypes+4);

using utils::RestrTraj;

class utils::RestrTraj_i: public gio::Ginstream
{
  friend class utils::RestrTraj;

  std::string d_current;
  //int d_switch;
  int d_skip;
  int d_stride;

  // the TIMESTEP block information
  bool d_time_read;
  int d_step;
  double d_time;

  // the JVALUERESEPS block information
  bool d_jvaluereseps_read;
  JValueRestrData m_jvalueresepsdata;

  //  the XRAYRVALUE block information
  bool d_xrayrvalue_read;
  XrayRestrData m_xrayrestrdata;


  RestrTraj_i(const std::string &name, int skip, int stride)
  : d_current(),
  //d_switch(),
  d_skip(skip),
  d_stride(stride),
  d_time_read(false),
  d_step(0),
  d_time(0.0),
  d_jvaluereseps_read(false),
  d_xrayrvalue_read(false){
    open(name);
    getline(d_current);
  }
  ~RestrTraj_i(){}

  // method
  void readTimestep();
  void readJvalueResEps();
  void readXrayRvalue();
};

// Constructors
RestrTraj::RestrTraj(const std::string &name, int skip, int stride)
  : d_skip(skip),
    d_stride(stride),
    d_stride_eof(false)
{
  d_this=new RestrTraj_i(name, skip, stride);
}

RestrTraj::RestrTraj(int skip, int stride)
  : d_this(NULL),
    d_skip(skip),
    d_stride(stride),
    d_stride_eof(false)
{
}

RestrTraj::~RestrTraj(){
  if(d_this)delete d_this;
}

void RestrTraj::open(const std::string &name){
  if(d_this){
    // recover skip and stride
    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
  }

  d_this=new RestrTraj_i(name, d_skip, d_stride);
}

void RestrTraj::close(){
  if(d_this){
    d_this->close();

    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
    d_this=NULL;
  }
}

bool RestrTraj::eof()const{
  if (!d_this){
    throw RestrTraj::Exception("eof, but no open file");
  }

  return d_this->stream().eof();
}

std::string RestrTraj::title()const{
  if (!d_this){
    throw RestrTraj::Exception("no open file");
  }
  return d_this->title();
}

void utils::RestrTraj_i::readTimestep()
{
  std::vector<std::string> buffer;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0) {
    throw RestrTraj::Exception("Coordinate file " + name() +
            " is corrupted. No END in TIMESTEP"
            " block. Got\n"
            + buffer[buffer.size() - 1]);
  }
  d_time_read = false;
  std::vector<std::string>::iterator it=buffer.begin();
  std::istringstream line(*it);

  line >> d_step >> d_time;
  if (line.fail()) {
    throw RestrTraj::Exception("Coordinate file " + name() +
            " is corrupted. Bad line in TIMESTEP block. Got\n" +
            *it);
  }
  // we read the time from the system
  d_time_read = true;
}

void utils::RestrTraj_i::readJvalueResEps(){
  // std::cerr << "readjvaluereseps" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it, to;
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw RestrTraj::Exception("Coordinate file " + name() +
			   " is corrupted. No END in JVALUERESEPS"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);

  d_jvaluereseps_read = false;
  m_jvalueresepsdata.data().clear();

  it=buffer.begin();
  to = buffer.end() - 1;

  for (; it != to; ++it) {
    _lineStream.clear();
    _lineStream.str(*it);

    JValueRestrData::JValueEps dat;
    _lineStream >> dat.i >> dat.j >> dat.k >> dat.l;
    if (_lineStream.fail())
      throw RestrTraj::Exception("Bad line in JVALUERESEPS block:\n" + *it +
            "\nTrying to read atom numbers");
    ++it;
    if (it == to)
      throw RestrTraj::Exception("Not enough lines in JVALUERESEPS block\n");
    _lineStream.clear();
    _lineStream.str(*it);

    double eps;
    while((_lineStream >> eps)) {
      dat.epsilon.push_back(eps);
    }

    m_jvalueresepsdata.data().push_back(dat);
  }
  d_jvaluereseps_read = true;
}

RestrTraj &RestrTraj::operator>>(JValueRestrData &data) {
  if (d_this->d_jvaluereseps_read)
    data = d_this->m_jvalueresepsdata;
  else
    throw RestrTraj::Exception("no jvalue epsilon data in traj file.");

  return *this;
}

void utils::RestrTraj_i::readXrayRvalue(){
  std::vector<std::string> buffer;
  getblock(buffer);

  d_xrayrvalue_read = false;

  if(buffer.size() != 9 || buffer[buffer.size()-1].find("END")!=0)
    throw RestrTraj::Exception("Coordinate file " + name() +
			   " XRAYRVALUE block is corrupted.");

  _lineStream.clear();
  _lineStream.str(buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4]
   + buffer[5] + buffer[6] + buffer[7]);

  _lineStream >> m_xrayrestrdata.state().scale_inst
          >> m_xrayrestrdata.state().r_inst
          >> m_xrayrestrdata.state().scale_free_inst
          >> m_xrayrestrdata.state().r_free_inst
          >> m_xrayrestrdata.state().scale_avg
          >> m_xrayrestrdata.state().r_avg
          >> m_xrayrestrdata.state().scale_free_avg
          >> m_xrayrestrdata.state().r_free_avg;

  if(_lineStream.fail())
    throw RestrTraj::Exception("Coordinate file " + name() +
			   ": bad line in XRAYRVALUE block.");

  d_xrayrvalue_read = true;
}

RestrTraj &RestrTraj::operator>>(XrayRestrData &data) {
  if (d_this->d_xrayrvalue_read)
    data = d_this->m_xrayrestrdata;
  else
    throw RestrTraj::Exception("no xray restraint data in traj file.");

  return *this;
}


void RestrTraj::read(){
  if (!d_this){
    throw RestrTraj::Exception("read in frame, but no open file");
  }

  if(!d_this->stream())
    throw Exception("File "+name()+" is corrupted.");

  const std::string first =d_this->d_current;
  // std::cerr << first << std::endl;
  std::vector<std::string> buffer;

  // skip frames
  // std::cerr << "operator<< : skip=" << d_this->d_skip << std::endl;

  for( ; d_this->d_skip > 0; --d_this->d_skip){

    do{
      // std::cerr << "skipping block " << d_this->d_current << std::endl;

      d_this->skipblock();
      d_this->getline(d_this->d_current);

    } while (d_this->d_current != first &&
	     (!d_this->stream().eof()));

    if (d_this->stream().eof()){
      // std::cerr << "skip eof: " << d_this->d_skip << std::endl;

      --d_this->d_skip;
      d_stride_eof = true;
      return;
    }

  }

  // only stride if not skip because of eof during last stride
  if (d_stride_eof == false){
    int i=1;
    for( ; i < d_this->d_stride; ++i){

      do{

	d_this->skipblock();
	d_this->getline(d_this->d_current);

      } while (d_this->d_current != first &&
	       (!d_this->stream().eof()));

      if (d_this->stream().eof()){
	// safe remaining strides in skip for next file

	std::cerr << "stride eof: " << d_this->d_stride
		  << "\ti: " << i << std::endl;

	d_this->d_skip = d_this->d_stride - i - 1;
	d_stride_eof = true;
	return;
      }

    }
  }

  // skipping and striding worked...
  d_stride_eof = false;

  do{
    switch(BLOCKTYPE[d_this->d_current]){
      case titleblock :
        break;
      case timestep:
	d_this->readTimestep();
	break;
      case jvaluereseps:
	d_this->readJvalueResEps();
	break;
      case xrayrvalue:
        d_this->readXrayRvalue();
        break;
      default:
	throw
	  Exception("Block "+d_this->d_current+
		    " is unknown in a restraint trajectory file");
	break;
    }
    d_this->getline(d_this->d_current);
  } while(d_this->d_current!=first&&!d_this->stream().eof());
  return;
}

RestrTraj &RestrTraj::operator>>(utils::Time &time){
  if (time.read()) {
    // we read the time from the traj
    if (!d_this->d_time_read) {
      throw Exception("trajectory does not contain a TIMESTEP block. "
              "Use @time to specify the time.");
    }
    time.time() = d_this->d_time;
  } else {
    // we have to calculate the time
    time.time() += time.dt();
  }

  return *this;
}

std::string RestrTraj::name()const{
  return d_this->name();
}

int RestrTraj::skip()const
{
  if (d_this)
    return d_this->d_skip;
  else return d_skip;
}

int RestrTraj::stride()const
{
  if (d_this)
    return d_this->d_stride;
  else return d_stride;
}

void RestrTraj::skip(int skip)
{
  if (d_this)
    d_this->d_skip = skip;
  else
    d_skip = skip;
}

void RestrTraj::stride(int stride)
{
  if (d_this)
    d_this->d_stride = stride;
  else
    d_stride = stride;
}

bool RestrTraj::stride_eof()const
{
  return d_stride_eof;
}
