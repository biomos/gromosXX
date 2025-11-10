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
#include "../gcore/BbLink.h"

namespace gcore{
  class BuildingBlock;
}

namespace gio{

  class InLinkingBlock_i;
  /**
   * Class InLinkingBlock
   * Instream that can read in a topology linking file
   *
   * @class InLinkingBlock
   * @author C. Oostenbrink
   * @ingroup gio
   * @sa gcore::BuildingBlock
   * @sa gcore::BbSolute
   + @sa gcore::BbLink
   * @sa gcore::SolventTopology
   */
  class InLinkingBlock{
    InLinkingBlock_i *d_this;
    // prevent default construction, copying and assignment
    InLinkingBlock();
    InLinkingBlock(const InLinkingBlock &);
    InLinkingBlock &operator=(const InLinkingBlock &);
    
  public:
    /**
     * open a building block file
     * @param file the file
     */
    InLinkingBlock(std::string file);
    
    ~InLinkingBlock();

    /**
     * access the building blocks
     */
    const std::vector<gcore::BbLink> &links()const;


    /**
     * access the force field code
     */
    const std::string forceField()const;

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
	gromos::Exception("InLinkingBlock", what){}
    };
  };
}
#endif
