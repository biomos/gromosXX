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

/* 	$Id$	 */

#ifndef MAXFLOAT
#define	MAXFLOAT	((float)3.40282346638528860e+38)
#endif

#ifndef INCLUDED_UTILS_PROPERTYCONTAINER
#define INCLUDED_UTILS_PROPERTYCONTAINER

#include <vector>
#include <string>
#include <iostream>

#include "Value.h"
#include "Property.h"
#include "../gmath/Stat.h"
#include "../gromos/Exception.h"
#include "../gmath/Distribution.h"

namespace gcore
{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace gmath
{
  class Distribution;
}

namespace bound
{
  class Boundary;
}

namespace utils
{
  /**
   * Class PropertyContainer
   * Purpose: Implements a container for properties.
   *
   * Description:
   * This class implements a container for properties. It parses the
   * arguments and constructs the single properties. It loops over these
   * properties to calculate their values, averages or to print out the
   * information.
   *
   * @class PropertyContainer
   * @version Wed Jul 31 2002
   * @author M. Christen
   * @ingroup utils
   * @sa utils::Property
   */
  
  class PropertyContainer: public std::vector<Property *>
  {
  public:
    /**
     * Constructor.
     */
    PropertyContainer();
    /**
     * Constructor.
     * Most properties need a reference to the system.
     */
    PropertyContainer(gcore::System &sys, bound::Boundary *pbc);
    /**
     * Destructor.
     */
    virtual ~PropertyContainer();
    
    /**
     * Clear the property container and assign new system and boundary
     */
    void reinitialize(gcore::System &sys, bound::Boundary *pbc);
    /**
     * Add a property specifier (thereby constructing a property class).
     */
    int addSpecifier(std::string s);
    /**
     * Get title string.
     */
    std::string toTitle()const;
    /**
     * Get results of a calculation.
     */
    std::string toString()const;
    /**
     * Calculate all properties in the container.
     */
    void calc();

    /**
     * @struct Exception
     * PropertyContainer exception.
     */
    struct Exception: public gromos::Exception
    {
      /**
       * Constructor.
       */
      Exception(const std::string &what) : 
	gromos::Exception("PropertyContainer", what) {}
    };

  protected:
    // internal functions
    /**
     * Parse the arguments string into property type and the rest of
     * the arguments.
     */
    void parse(std::string s);
    /**
     * separate s by ';' and parse the properties
     */
    void parse_multiple(std::string s, std::vector<Property *> & prop);
    /**
     * Parse the arguments string into property type and parse the rest of
     * the properties arguments.
     */
    void parse_single(std::string s, std::vector<Property *> & prop);
    /**
     * create an AverageProperty and add the properties specified in
     * the angular (<>) brackets.
     */
    void parse_average(std::string s);
    /**
     * Creates the specified property. Override if you want to add user
     * defined properties (ie properties derived from Property, but not
     * DistanceProperty, AngleProperty or TorsionProperty).
     * The remaining arguments are passed to the property to parse.
     */      
    virtual Property *createProperty
    (
     std::string type, 
     std::vector<std::string> const & arguments,
     int x
     );

    /**
     * Reference to the system. Often needed to construct a property.
     */
    gcore::System *d_sys;
    /**
     * and the boundary conditions
     */
    bound::Boundary *d_pbc;

    /**
     * Prints all the values of the properties in the container.
     */
    friend std::ostream & operator<<(std::ostream &os, PropertyContainer const & s);
  };  
}

#endif

  
