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

// utils_DipTraj.cc
#include "DipTraj.h"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "../gmath/Vec.h"
#include "../gio/Ginstream.h"
#include "groTime.h"

enum blocktype {
  titleblock, timestep, dip
};

typedef std::map<std::string, blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TITLE", titleblock),
  BT("TIMESTEP", timestep),
  BT("DIPOLE", dip)};


static std::map<std::string, blocktype> BLOCKTYPE(blocktypes, blocktypes + 3);

using utils::DipTraj;

class utils::DipTraj_i : public gio::Ginstream {
  friend class utils::DipTraj;

  std::string d_current;
  //int d_switch;
  int d_skip;
  int d_stride;

  // the TIMESTEP block information
  bool d_time_read;
  int d_step;
  double d_time;

  // the DIPOLE block information
  bool d_dip_read;
  utils::DipData m_dipdata;

  DipTraj_i(const std::string &name, int skip, int stride)
  : d_current(),
  //d_switch(),
  d_skip(skip),
  d_stride(stride),
  d_time_read(false),
  d_step(0),
  d_time(0.0),
  d_dip_read(false) {
    open(name);
    getline(d_current);
  }

  ~DipTraj_i() {
  }

  // method
  void readTimestep();
  void readDipole();
};

// Constructors

DipTraj::DipTraj(const std::string &name, int skip, int stride)
: d_skip(skip),
d_stride(stride),
d_stride_eof(false) {
  d_this = new DipTraj_i(name, skip, stride);
}

DipTraj::DipTraj(int skip, int stride)
: d_this(NULL),
d_skip(skip),
d_stride(stride),
d_stride_eof(false) {
}

DipTraj::~DipTraj() {
  if (d_this)delete d_this;
}

void DipTraj::open(const std::string &name) {
  if (d_this) {
    // recover skip and stride
    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
  }

  d_this = new DipTraj_i(name, d_skip, d_stride);
}

void DipTraj::close() {
  if (d_this) {
    d_this->close();

    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
    d_this = NULL;
  }
}

bool DipTraj::eof()const {
  if (!d_this) {
    throw DipTraj::Exception("eof, but no open file");
  }

  return d_this->stream().eof();
}

std::string DipTraj::title()const {
  if (!d_this) {
    throw DipTraj::Exception("no open file");
  }
  return d_this->title();
}

void utils::DipTraj_i::readTimestep() {
  std::vector<std::string> buffer;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0) {
    throw DipTraj::Exception("Coordinate file " + name() +
            " is corrupted. No END in TIMESTEP"
            " block. Got\n"
            + buffer[buffer.size() - 1]);
  }
  d_time_read = false;
  std::vector<std::string>::iterator it = buffer.begin();
  std::istringstream line(*it);

  line >> d_step >> d_time;
  if (line.fail()) {
    throw DipTraj::Exception("Coordinate file " + name() +
            " is corrupted. Bad line in TIMESTEP block. Got\n" +
            *it);
  }
  // we read the time from the system
  d_time_read = true;
}

void utils::DipTraj_i::readDipole() {
  //std::cerr << "readdip" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it, to;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw DipTraj::Exception("Coordinate file " + name() +
          " is corrupted. No END in DIPOLE"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  d_dip_read = false;
  m_dipdata.data().clear();

  it = buffer.begin();
  to = buffer.end() - 1;

  for (; it != to; ++it) {
    _lineStream.clear();
    _lineStream.str(*it);

    utils::DipData::Dip dat;
    _lineStream >> dat.dipole[0] >> dat.dipole[1] >> dat.dipole[2] >> dat.vol;
    if (_lineStream.fail())
      throw DipTraj::Exception("Bad line in DIPOLE block:\n" + *it +
            "\nTrying to read atom numbers");

    m_dipdata.data().push_back(dat);
  }
  d_dip_read = true;
}

DipTraj &DipTraj::operator>>(DipData &data) {
  if (d_this->d_dip_read)
    data = d_this->m_dipdata;
  else
    throw DipTraj::Exception("no dipole moment data in traj file.");

  return *this;
}

void DipTraj::read() {
  if (!d_this) {
    throw DipTraj::Exception("read in frame, but no open file");
  }

  if (!d_this->stream())
    throw Exception("File " + name() + " is corrupted.");

  const std::string first = d_this->d_current;
  //std::cerr << first << std::endl;
  std::vector<std::string> buffer;

  // skip frames
  // std::cerr << "operator<< : skip=" << d_this->d_skip << std::endl;

  for (; d_this->d_skip > 0; --d_this->d_skip) {
    do {
      // std::cerr << "skipping block " << d_this->d_current << std::endl;
      d_this->skipblock();
      d_this->getline(d_this->d_current);
    } while (d_this->d_current != first && (!d_this->stream().eof()));

    if (d_this->stream().eof()) {
      // std::cerr << "skip eof: " << d_this->d_skip << std::endl;
      --d_this->d_skip;
      d_stride_eof = true;
      return;
    }

  }

  // only stride if not skip because of eof during last stride
  if (d_stride_eof == false) {
    int i = 1;
    for (; i < d_this->d_stride; ++i) {
      do {
        d_this->skipblock();
        d_this->getline(d_this->d_current);

      } while (d_this->d_current != first && (!d_this->stream().eof()));

      if (d_this->stream().eof()) {
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

  do {
    switch (BLOCKTYPE[d_this->d_current]) {
      case titleblock:
        break;
      case timestep:
        d_this->readTimestep();
        break;
      case dip:
        d_this->readDipole();
        break;
      default:
        throw
        Exception("Block " + d_this->d_current +
                " is unknown in a dipole trajectory file");
        break;
    }
    d_this->getline(d_this->d_current);
  } while (d_this->d_current != first && (!d_this->stream().eof()));
  return;
}

DipTraj &DipTraj::operator>>(utils::Time &time) {
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

std::string DipTraj::name()const {
  return d_this->name();
}

int DipTraj::skip()const {
  if (d_this)
    return d_this->d_skip;
  else return d_skip;
}

int DipTraj::stride()const {
  if (d_this)
    return d_this->d_stride;
  else return d_stride;
}

void DipTraj::skip(int skip) {
  if (d_this)
    d_this->d_skip = skip;
  else
    d_skip = skip;
}

void DipTraj::stride(int stride) {
  if (d_this)
    d_this->d_stride = stride;
  else
    d_stride = stride;
}

bool DipTraj::stride_eof()const {
  return d_stride_eof;
}
