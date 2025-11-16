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
 * @file InCIF.h
 * reading a cif file
 */

#ifndef INCLUDED_INCIF_H
#define	INCLUDED_INCIF_H

#include <string>
#include <vector>
#include <fstream>

namespace gio {

  /**
   * @struct CIFData
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief container for CIF data
   *
   * Data type to hold a structure factor read from a crystallographic
   * information file (CIF).
   */
  struct CIFData {
    /**
     * h Miller index
     */
    int h;
    /**
     * k Miller index
     */
    int k;
    /**
     * l Miller index
     */
    int l;
    /**
     * observed structure factor amplitude
     */
    double f_obs;
    /**
     * standard deviation of observed structure factor amplitude
     */
    double stddev_f_obs;
  };

  /**
   * @class InCIF
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief reads a CIF file
   * A class to handle and read crystallographic information files (CIF).
   */
  class InCIF {
  public:
    /**
     * constructor
     */
    InCIF() {}
    /**
     * construct from file name
     * @param file file name
     */
    InCIF(std::string file);
    /**
     * destructor
     */
    ~InCIF();
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
     * check whether file is open
     * @return true if the file is open, false otherwise.
     */
    bool isOpen();
    /**
     * get the cif data
     * @return a vector containing the CIF data
     */
    std::vector<CIFData> getData();
  protected:
    std::ifstream file_stream;
  };
}

#endif	/* _INCIF_H */

