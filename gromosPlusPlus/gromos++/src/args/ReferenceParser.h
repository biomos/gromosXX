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

// args_ReferenceParser.h

#ifndef INCLUDED_ARGS_REFERENCEPARSER
#define INCLUDED_ARGS_REFERENCEPARSER

#include <vector>
#include "../fit/Reference.h"

namespace gcore{
  class System;
}

namespace args{
  class Arguments;
}

namespace fit{
  class Reference;
}

namespace args{
  /**
   * Class ReferenceParser
   * Parses the different inputs for specifying the reference from the 
   * arguments 
   * 
   * The reference parser looks at a set of input flags that can be used
   * to specify the reference atoms (class, atoms). These are required 
   * for a (rotational) fit. The atoms are added to the reference
   * immediately.
   *
   * @class ReferenceParser
   * @author V. Kraeutler
   * @ingroup args
   * @sa args::Arguments
   */
  class ReferenceParser{
    
    std::vector<int> myMols;
    gcore::System &mySys;
    const args::Arguments &myArgs;
    fit::Reference &myRef;
    bool myAdded;

  public:

    /**
     * ReferenceParser constructor
     *
     * @param sys The system on which the reference is based
     * @param args All the arguments, the parser only looks at the 
     *             the interesting ones
     * @param ref The reference to which the atoms need to be added
     */
    ReferenceParser(
      gcore::System &sys,
      const args::Arguments &args,
      fit::Reference &ref);
    /**
     * Function that is only needed internally
     */
    void add_ref();
    /**
     * Function that is only needed internally
     */
    void getMolecules();
    /**
     * Function that is only needed internally
     */
    void classes();
    /**
     * Function that is only needed internally
     */
    void atoms();
  };

}
#endif
