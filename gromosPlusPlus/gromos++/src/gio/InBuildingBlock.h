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

// gio_InBuildingBlock.h

#ifndef INCLUDED_GIO_IBUILDING
#define INCLUDED_GIO_IBUILDING

#include <string>

#include "../gromos/Exception.h"

namespace gcore{
  class BuildingBlock;
}

namespace gio{

  class InBuildingBlock_i;
  /**
   * Class InBuildingBlock
   * Instream that can read in a GROMOS mtb-file
   *
   * In addition to the standard blocks in the mtb-file, this stream can also
   * read in end-group blocks.
   * These blocks should appear at the end of the file, after the 
   * solvent blocks
   *
   * @class InBuildingBlock
   * @author B.C. Oostenbrink
   * @ingroup gio
   * @sa gcore::BuildingBlock
   * @sa gcore::BbSolute
   * @sa gcore::SolventTopology
   */
  class InBuildingBlock{
    InBuildingBlock_i *d_this;
    // prevent default construction, copying and assignment
    InBuildingBlock();
    InBuildingBlock(const InBuildingBlock &);
    InBuildingBlock &operator=(const InBuildingBlock &);
    
  public:
    /**
     * open a building block file
     * @param file the file
     */
    InBuildingBlock(std::string file);
    
    ~InBuildingBlock();

    /**
     * access the building blocks
     */
    const gcore::BuildingBlock &building()const;

    /**
     * access the title
     */
    const std::string title()const;

    /**
     * @struct Exception
     * @ingroup gio::InBuildingBlock
     * The exception type for BB input exceptions
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InBuildingBlock", what){}
    };
  };
}
#endif
