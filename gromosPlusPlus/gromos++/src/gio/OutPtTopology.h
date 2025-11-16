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

// gio_OutPtTopology.h

#ifndef INCLUDED_OUTPTTOPOLOGY
#define INCLUDED_OUTPTTOPOLOGY

#include<string>
#include <set>

namespace gcore{
  class PtTopology;
  class System;
}

namespace gio{
  /**
   * Class OutPtTopology
   * an out stream that defines how a perturbation topology should be written out
   *
   * @class OutPtTopology
   * @author N. Schmid
   * @ingroup gio
   */
  class OutPtTopology{
    std::string d_title;
    std::ostream &d_os;
  public:
    /**
     * construct the OutPtTopology from a stream
     * @param os the out stream to which the perturbation topology is written
     */
    OutPtTopology(std::ostream &os);
    /**
     * destruct the OutPtTopology
     */
    ~OutPtTopology();
    /**
     * sets the content of the title block
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * writes out a perturbation topology between two states A and B.
     * @param pt the perturbation topology
     * @param sys the reference system to detect the hydrogens. If NULL 
     *        (default) partitioning into hydrogen and non-hydrogen blocks
     *        is omitted.
     * @param a the index of state A (default 0)
     * @param b the index of state B (default 1)
     */
    void write(const gcore::PtTopology & pt,
               gcore::System * sys = NULL, int a=0, int b=1);
    /**
     * writes out a multiple perturbation topology (nonbonded only)
     * @param pt the perturbation topology
     */
    void write_multiple(const gcore::PtTopology & pt);
  };
}

#endif

