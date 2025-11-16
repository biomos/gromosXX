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
 * @file InIACElementNameMapping.h
 * a file to read element mapping files
 */

#ifndef INCLUDED_INIACELEMENTNAMEMAPPING_H
#define	INCLUDED_INIACELEMENTNAMEMAPPING_H

#include <string>
#include <map>

#include "Ginstream.h"

namespace gio {

  /**
   * @class InIACElementNameMapping
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief reads a IAC to element mapping file
   * A class to handle and read file that map IACs to element names
   *
   * Format of the IAC to element-name mapping file:
   * @verbatim
TITLE
map IAC to element name for xxxx force-field
END
ELEMENTMAPPING
# IAC ELEMENTNAME
1 O
11 N
12 C
26 Fe
END
@endverbatim
   */
  class InIACElementNameMapping {
  public:
    /**
     * constructor
     */
    InIACElementNameMapping() {}
    /**
     * construct from file name
     * @param file file name
     */
    InIACElementNameMapping(std::string file);
    /**
     * destructor
     */
    ~InIACElementNameMapping();
    /**
     * open a file
     * @param file file name
     */
    void open(std::string file);
    /**
     * close the file
     */
    void close();
    /**
     * get the mapping data
     * @return a map containing the IAC to element mapping
     */
    std::map<int, std::string> getData();
  protected:
    Ginstream file_stream;
  };
}


#endif	/* INCLUDED_INIACELEMENTNAMEMAPPING_H */

