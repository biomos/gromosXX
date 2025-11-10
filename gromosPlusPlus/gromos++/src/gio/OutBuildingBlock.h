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

// gio_OutBuildingBlock.h
#ifndef INCLUDED_OUTBUILDINGBLOCK
#define INCLUDED_OUTBUILDINGBLOCK

#include <string>
#include <set>
#include "../gcore/SolventTopology.h"

namespace gcore{
  class BuildingBlock;
  class BbSolute;
}

namespace gio{
  /**
   * Class OutBuildingBlock
   * an outstream that defines how a GROMOS molecular building block file should be written out
   *
   * @class OutBuildingBlock
   * @author N. Schmid
   * @ingroup gio
   */
  class OutBuildingBlock{
    std::string d_title;
    std::ostream &d_os;

    // prevent copying and assignment
    OutBuildingBlock();
    OutBuildingBlock(const OutBuildingBlock&);
    OutBuildingBlock &operator=(const OutBuildingBlock&);
  public:
    /**
     * The type of a building block
     */
    enum BBType { BBTypeSolute, BBTypeSolvent, BBTypeEnd };
    /**
     * create from a stream
     */
    OutBuildingBlock(std::ostream &os);
    ~OutBuildingBlock();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the BuildingBlock as a molecular topology building block file
     * @param bb the building block data
     */
    void write(const gcore::BuildingBlock &bb);
    /**
     * write a single building block 
     * @param bb the building block
     * @param type the type (either solute or end group)
     */
    void writeSingle(const gcore::BbSolute & bb, BBType type);
    /**
     * write a single solvent building block
     * @param bb the solvent building block
     */
    void writeSolvent(const gcore::SolventTopology & bb);
  };
}

#endif

