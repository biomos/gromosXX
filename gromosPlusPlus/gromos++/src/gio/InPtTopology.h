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

// gio_InTopology.h

#ifndef INCLUDED_GIO_INPTTOPOLOGY
#define INCLUDED_GIO_INPTTOPOLOGY

#include <string>

#include "../gromos/Exception.h"

namespace gcore{
  class System;
  class PtTopology;
}

namespace gio{
  class InPtTopology_i;
  /**
   * Class InPtTopology
   * defines an instream that can read in a perturbation topology
   *
   * @class InPtTopology
   * @ingroup gio
   * @author B.C. Oostenbrink  N. Schmid
   * @sa gcore::System
   * @sa gcore::PtTopology
   */
  class InPtTopology{
    InPtTopology_i *d_this;
    
  public:
    /**
     * Construct the InPtTopology from a file name and parse the content
     */
    InPtTopology(std::string str);
    /**
     * Destructor
     */
    ~InPtTopology();
    
    /**
     * Accessor to the perturbation topology read
     */
    const gcore::PtTopology &ptTopo()const;
    /**
     * Accessor to the title
     */
    const std::string title()const;

    /**
     * @struct Exception
     * The Exception that occurs when reading a perturbation topology
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InPtTopology", what){}
    };
  };
}
#endif
