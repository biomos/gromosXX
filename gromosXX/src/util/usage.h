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
 * @file usage.h
 */

#ifndef INCLUDED_USAGE_H
#define INCLUDED_USAGE_H

namespace util
{
  /**
   * @class Known
   * Known
   */
  class Known : public std::set<std::string>
  {
  public:
    Known & operator << (std::string s)
    {
      insert(s);
      return *this;
    }
  };

  void get_usage(Known const &knowns, std::string &s, std::string name);

  void print_title(bool mpi, std::ostream & os = std::cout,
          bool repex = false);
  
}

#endif

  
