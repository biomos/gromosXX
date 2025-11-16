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
 * @file groTime.cc
 * implements utils::Time
 */
#include "groTime.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "../args/Arguments.h"
#include "../gromos/Exception.h"

namespace utils {
  Time::Time(const args::Arguments & args) :
  d_t0(0.0), d_dt(1.0), d_current_time(0.0), d_steps(0), d_read(true), d_do_timeseries(false) {
    // get the time command line argument
    if (args.count("time") != -1)
      doSeries() = true;
    
    args::Arguments::const_iterator iter = args.lower_bound("time");
    if (iter != args.upper_bound("time")) {
      // argument supplied: do not read
      read() = false;
      std::istringstream in(iter->second);
      if (!(in >> start_time()))
        throw gromos::Exception("time", "@time: starting time is not numeric.");
      ++iter;
    }
    if (iter != args.upper_bound("time")) {
      std::istringstream in(iter->second);
      if (!(in >> dt()))
        throw gromos::Exception("time", "@time: timestep is not numeric.");
    }

    if (read() == false) {
      // the first frame has time start_time(). so when we call the operator
      // to read the time it will be increased to start_time() + dt(). That's
      // why we have to start at start_time() - dt()
      time() = start_time() - dt();
      steps() = 0.0; 
    }
  }

  void Time::print(std::ostream & out) const {
    // write the time
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setw(15) << std::setprecision(9) << time();
  }
  
  std::ostream & operator<<(std::ostream &os, Time const & t) {
    t.print(os);
    return os;
  }
}




