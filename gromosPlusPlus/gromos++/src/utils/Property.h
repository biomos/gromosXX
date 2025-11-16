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
 * @file Property.h
 * Properties and PropertySpecifiers
 */

#ifndef MAXFLOAT
#define	MAXFLOAT	((float)3.40282346638528860e+38)
#endif

#ifndef INCLUDED_UTILS_PROPERTY
#define INCLUDED_UTILS_PROPERTY

#include <vector>
#include <string>

#include "Value.h"
#include "VectorSpecifier.h"
#include "ExpressionParser.h"
#include "../gmath/Stat.h"
#include "../gromos/Exception.h"

namespace gcore
{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace bound
{
  class Boundary;
}

namespace utils
{

  /**
   * Class Property 
   * Purpose: General base class for properties to be calculated
   * during an analysis run.
   *
   * Description:
   * This class defines (and implements) the general methods that a
   * property calculated during an analysis run should posses.
   *
   * @section PropertySpecifier Property Specifier
   * A property specifier is of the general form:<br>
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <type>%<arg1>[%<arg2>...] @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">type</span> can be:
   *   - <b>d</b> @ref DistanceProperty "Distance"
   *   - <b>a</b> @ref AngleProperty "Angle"
   *   - <b>t</b> @ref TorsionProperty "Torsion"
   *   - <b>tp</b> @ref PeriodicTorsionProperty "periodic Torsion, as TorsionProperty but mapped to -180 to 180 degrees"
   *   - <b>ct</b> @ref CrossTorsionProperty "CrossTorsion"
   *   - <b>hb</b> @ref HBProperty "Hydrogen bond"
   *   - <b>st</b> @ref StackingProperty "Stacking"
   *   - <b>o</b> @ref VectorOrderProperty "Order"
   *   - <b>op</b> @ref VectorOrderParamProperty "Order parameter"
   *   - <b>pr</b> @ref PseudoRotationProperty "Pseudo rotation"
   *   - <b>pa</b> @ref PuckerAmplitudeProperty "Pucker amplitude"
   *   - <b>expr</b> @ref ExpressionProperty "Expression Property"
   * - <span style="color:darkred;font-family:monospace">&lt;arg1&gt;</span> is usualy
   *   and @ref AtomSpecifier or @ref VectorSpecifier
   * - further arguments are depdendent on the type.
   *
   * @subsection property_av_dist Averages and Distributions
   *
   * It is also possible to calculate the @ref AverageProperty "average and 
   * distribution" of a set of properties.
   *
   * @subsection property_subst Variable Substitution
   * If 'a' is specified for the molecules, a property for each of the molecules is generated.
   * If 'a' is specified for an atom, all atoms of the molecule are taken.
   * Instead of 'a' also 'x' may be used. The values inserted for 'x' are specified as follows:
   * @verbatim 
   @prop 'd%(x):1,2|x=345,350-360'
   @prop 'd%1:(x),(x+1)|x=3,8'
   @endverbatim
   *
   * The first will create distance properties (between atom 1 and 2) of molecules 345 and 350 to 360.
   * The second generates two properties, the first is d%1:3,4, the second d%1:8,9.
   * Note that the brackets are necessary around expressions (to distinguish a '-' sign from a range).
   * 
   * @class Property
   * @version Jul 31 15:00 2002
   * @author M. Christen
   * @ingroup utils
   * @sa AtomSpecifier VectorSpecifier utils::PropertyContainer utils::DistanceProperty
   *     utils::AngleProperty utils::TorsionProperty utils::CrossTorsionProperty utils::HBProperty
   *     utils::StackingProperty utils::VectorOrderProperty
   *     utils::VectorOrderParamProperty utils::PseudoRotationProperty 
   *     utils::PuckerAmplitudeProperty utils::ExpressionProperty
   */
  class Property;
  std::ostream &operator<<(std::ostream &os, Property const & p);
  class Property
  {
  public:
    /**
     * Maximum number of arguments.
     */
    static const int MAXARGUMENTS = 10;

    /**
     * Standard constructor.
     * Takes a reference to the system, so that the properties can
     * calculate their value (ie by looking up the atom positions).
     */
    Property(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~Property();

    // accessor
    std::string type() const {return d_type;}

    /**
     * Return the value of the property.
     */
    Value const & getValue() const { return d_value;
    }

    /**
     * Return the i-th value of the property.
     */
    Value getValue(int i) const {
      switch (d_value.type()) {
        case val_scalar:
          return Value(getScalarStat().data()[i]);
        case val_vecspec:
        case val_vector:
          return Value(getVectorStat().data()[i]);
      }
      return Value(0.0);
    }
    /**
     * return the number of values contained
     * @return the number of values
     */
    int num() const {
      switch(d_value.type()) {
        case val_scalar :
          return getScalarStat().n();
        case val_vecspec:
        case val_vector:
          return getVectorStat().n();
      }
      return 0;
    }


    /**
     * arguments accessor
     */
    std::vector<Value> const & args() const { return d_arg; }
    
    /**
     * As most of the properties i can think of have something to do with
     * atoms and molecules, i define these vectors in the base class.
     * This is also done in order to be able to write one general
     * arguments parsing function in the base class.
     */
    AtomSpecifier & atoms() { return d_atom; }
    /**
     * and a const variant
     */
    AtomSpecifier const & atoms()const { return d_atom; }
    /**
     * scalar stat
     */
    const gmath::Stat<double> & getScalarStat() const { return d_scalar_stat; }
    /**
     * vector stat
     */
    const gmath::Stat<gmath::Vec> & getVectorStat() const { return d_vector_stat; }
    /**
     * scalar stat
     */
    gmath::Stat<double> & getScalarStat() { return d_scalar_stat; }
    /**
     * vector stat
     */
    gmath::Stat<gmath::Vec> & getVectorStat() { return d_vector_stat; }
    
    // methods
    
    /**
     * Calculates the value of the property.
     * Override in derived classes.
     */
    virtual Value const & calc() = 0;
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;
    /**
     * Returns the value in string representation.
     */
    virtual std::string toString()const { return d_value.toString(); }
    /**
     * return the average, rmsd and error estimate of all calculations so far
     */
    virtual std::string average()const;    
    /**
     * Returns the type of the interaction from the
     * topology.
     */
    virtual int getTopologyType(gcore::System const &sys);
    /**
     * @struct Exception
     * Property exception
     */
    struct Exception: public gromos::Exception
    {
      /**
       * Constructor.
       */
      Exception(const std::string &what): 
	gromos::Exception("Property", what) {}
    };
    /**
     * Parse the command line property specification.
     * This is the standard implementation. It knows about
     * atoms and arguments
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);

    /**
     * nearest image distance between to values. i.e. a "vector" pointing from
     * the first to the second value.
     * @param first value
     * @param second value
     * @return the nearest image distance value
     */
    virtual Value nearestImageDistance(const Value & first, const Value & second) const;

  protected:
    /**
     * Helper method to parse the atom part of the property specifier.
     * The argument is in AtomSpecifier format.
     */
    void parseAtoms(std::string atoms, int x);
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * add a value to the stats
     */
    void addValue(Value const & v);
    
    // member variables
    /**
     * Number of required arguments. Used in the standard parse function.
     * Set in the constructor.
     */
    int REQUIREDARGUMENTS; 
    /**
     * The property type (in string representation).
     */
    std::string d_type;
    /**
     * The atoms belonging to this property.
     */
    AtomSpecifier d_atom;
    /**
     * the current value
     */
    Value d_value;
    /**
     * the ideal (zero) value
     */
    std::vector<Value> d_arg;
    /**
     * Reference of the system.
     */
    gcore::System *d_sys;
    /**
     * pointer to a boundary object
     */
    bound::Boundary *d_pbc;

    /**
     * Statistics of all calculated scalar values
     */
    gmath::Stat<double> d_scalar_stat;

    /**
     * Statistics of all calculated vector values
     */
    gmath::Stat<gmath::Vec> d_vector_stat;

  };


  ////////////////////////////////////////////////////////////////////////////////

  /**
   * Class AverageProperty
   * Purpose: Meta property that averages over a set of properties
   *
   * Description:
   * The AverageProperty class provides a 'meta-' property, a property
   * that contains a list of other properties and averages over those.
   *
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <<prop1> <prop2>[ <props>...]> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<prop1\></span> and
   *   <span style="color:darkred;font-family:monospace">\<prop2\></span> are a
   *   @ref Property
   *
   * This calculates the average, rmsd and error estimate of the properties
   * <span style="color:darkred;font-family:monospace">\<prop1\></span>,
   * <span style="color:darkred;font-family:monospace">\<prop2\></span>... given 
   * between the angular brackets ('<', '>'). Multiple properties are seperated
   * by a space (' ').
   *
   * For example:
   * - @verbatim <d%1:1,10 d%2:1,10> @endverbatim means the average of the 
   *  distances between atoms 1 and 10 of the first and second molecule.
   *
   * Subtitution can be used within the averaging brackets.
   *
   * For example:
   * - @verbatim <d%x:1,10|x=1,2> @endverbatim calculates the same average
   *   as the example above.
   *
   * @class AverageProperty
   * @ingroup utils
   * @version Tue Aug 23 2005
   * @author markus
   * @sa utils::Property gmath::Stat
   */

  class AverageProperty : public Property
  {    
  public:
    /**
     * Constructor.
     */
    AverageProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~AverageProperty() {}
    /**
     * add a property
     */
    void addProperty(Property *p)
    {
      d_property.push_back(p);
    }
    /**
     * property accessor
     */
    Property * property(unsigned int i)
    {
      assert(i < d_property.size());
      return d_property[i];
    }
    /**
     * const property accessor
     */
    const Property * property(unsigned int i)const
    {
      assert(i < d_property.size());
      return d_property[i];
    }
    /**
     * properties accessor
     */
    std::vector<Property *> & properties()
    {
      return d_property;
    }
    /**
     * Calculate all properties.
     */
    virtual Value const & calc();
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;
    /**
     * Returns the value in string representation.
     */
    virtual std::string toString()const;

  protected:
    /**
     * the properties to average over
     */
    std::vector<Property *> d_property;
    /**
     * single calc scalar statistics
     */
    gmath::Stat<double> d_single_scalar_stat;
    /**
     * single calc vector statistics
     */
    gmath::Stat<gmath::Vec> d_single_vector_stat;

  };

  ////////////////////////////////////////////////////////////////////////////////

  /**
   * Class DistanceProperty
   * Purpose: Implements a distance property.
   *
   * Description:
   * This class implements a distance property. It is derived from the 
   * Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim d%<atomspec>%<step>%<lower_bound>%<upper_bound> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">\<step\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<lower_bound\></span> and
   *   <span style="color:darkred;font-family:monospace">\<upper_bound\></span> 
   *   are additional parameters
   *
   *
   * This calculates the distances between the two atoms in the atom specifier 
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span>. The
   * periodic boundary conditions are taken into account.
   *
   * The additional parameters are only used by certain @ref programs "programs"
   * like @ref gca. Their meaning may change depdending on the program.
   *
   * For example:
   * - @verbatim d%1:1,2 @endverbatim means the distance between atoms 1 and 2
   *   of molecule 1.
   * - @verbatim d%1:1;2:1 @endverbatim means the distance between the first 
   *   atoms of molecule 1 and 2.
   * - @verbatim d%va(com,1:a);va(com,2:a) @endverbatim means the distance
   *   between the centres of mass of molecules 1 and 2.
   *
   * Multiple distances are available via substitution.
   * 
   * For example:
   * - @verbatim d%1:res((x):H,N)|x=3-5 @endverbatim means the distances between
   *   the H and N atoms of residues 3 to 5 of molecule 1.
   *
   * @class DistanceProperty
   * @ingroup utils
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class DistanceProperty : public Property
  {    
  public:
    /**
     * Constructor.
     */
    DistanceProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~DistanceProperty() {}
    
    /**
     * Parses the arguments. Is overridden to check the input, but basically
     * just calls Property::parse.
     * @sa utils::Property::parse
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);
    /**
     * Calculates the distance.
     */
    virtual Value const & calc();

  protected:

    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
  };

  /**
   * Class AngleProperty
   * Purpose: Implements an angle property.
   *
   * Description:
   * Implements an angle property. It is derived from the Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim a%<atomspec>%<step>%<lower_bound>%<upper_bound> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">\<step\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<lower_bound\></span> and
   *   <span style="color:darkred;font-family:monospace">\<upper_bound\></span> 
   *   are additional parameters
   *
   * This calculates the angle between the three atoms defined by the
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span>
   * atom specifier. The periodic boundary conditions are taken into account.
   *
   * The additional parameters are only used by certain @ref programs "programs"
   * like @ref gca. Their meaning may change depdending on the program.
   *
   * For example:
   * - @verbatim a%1:1-3 @endverbatim means the angle defined by atoms 1, 2 and
   *   3 of molecule 1.
   * - @verbatim a%1:1,2;2:1 @endverbatim  means the angle defined by atom 1 and
   *   2 of the first molecule and atom 1 of the second molecule.
   * - @verbatim a%va(1,1:res(1:N,CA,CB,C));1:res(1:CA,C) @endverbatim means the
   *   angle between the virtual alpha hydrogen of residue 1 and the atoms CA 
   *   and C of residue 1 of molecule 1.
   *
   * Multiple angles are available via substitution.
   * 
   * For example:
   * - @verbatim a%1:(x),(x+1),(x+2)|x=1,2 @endverbatim means the angles defined
   *   by atoms 1-3 and 2-4 of molecule 1.
   *
   * @class AngleProperty
   * @ingroup utils
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class AngleProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    AngleProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~AngleProperty() {}
    
    /**
     * Parse and check property specifier (given in arguments).
     * Calls Property::parse and checks the arguments.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);
    /**
     * Calculate the angle between the given atoms.
     */
    virtual Value const & calc();
    /**
     * calculate the nearest image distance
     */
    virtual Value nearestImageDistance(const Value & first, const Value & second) const;
  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
  };
      
  /**
   * Class TorsionProperty
   * Purpose: Implements a torsion property.
   *
   * Description:
   * This class implements a torsion property. It is derived from the 
   * Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim t%<atomspec>%<step>%<lower_bound>%<upper_bound> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">\<step\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<lower_bound\></span> and
   *   <span style="color:darkred;font-family:monospace">\<upper_bound\></span> 
   *   are additional parameters
   *
   * This calculates the torsion definied by the four atoms in the 
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span>
   * atom specifier. The periodic boundary conditions are taken into account.
   *
   * The additional parameters are only used by certain @ref programs "programs"
   * like @ref gca. Their meaning may change depdending on the program.
   *
   * For example:
   * - @verbatim t%1:2-5 @endverbatim means the torsion defined by the atoms 2,
   *   3, 4 and 5 of molecule 1.
   * - @verbatim t%1:res(2:H,N,CA,CB) @endverbatim means the torsion defined by
   *   the atoms H, N, CA and CB of residue 2 of molecule 1.
   *
   * Multiple torsions are available via substitution.
   *
   * For example:
   * - @verbatim t%1:res((x):H,N,CA,C)|x=3-5 @endverbatim means the torsions 
   * defined bt the atoms H, N, CA and C of residues 3 to 5 of molecule 1.
   *
   * @class TorsionProperty
   * @ingroup utils
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */

  class TorsionProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    TorsionProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~TorsionProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);
    /**
     * Calculate the torsional angle.
     */
    virtual Value const & calc();
    /**
     * calculate the nearest image distance
     */
    virtual Value nearestImageDistance(const Value & first, const Value & second) const;
  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
  };
  
  /**
   * Class PeriodicTorsionProperty
   * Purpose: Implements a torsion property.
   *
   * Description:
   * This class implements a periodic torsion property. It is derived from the 
   * Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim t%<atomspec>%<step>%<lower_bound>%<upper_bound> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">\<step\></span>, 
   *   <span style="color:darkred;font-family:monospace">\<lower_bound\></span> and
   *   <span style="color:darkred;font-family:monospace">\<upper_bound\></span> 
   *   are additional parameters
   *
   * This calculates the periodic torsion definied by the four atoms in the 
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span>
   * atom specifier. The periodic boundary conditions are taken into account.
   *
   * The additional parameters are only used by certain @ref programs "programs"
   * like @ref gca. Their meaning may change depending on the program.
   *
   * For example:
   * - @verbatim t%1:2-5 @endverbatim means the torsion defined by the atoms 2,
   *   3, 4 and 5 of molecule 1.
   * - @verbatim t%1:res(2:H,N,CA,CB) @endverbatim means the torsion defined by
   *   the atoms H, N, CA and CB of residue 2 of molecule 1.
   *
   * Multiple torsions are available via substitution.
   * 
   * This class ist different to @ref TorsionProperty in the way taht the dihedral
   * angles are not summed up to follow the real rotations, but are mapped into
   * the range of -pi to pi.
   *
   * For example:
   * - @verbatim t%1:res((x):H,N,CA,C)|x=3-5 @endverbatim means the torsions 
   * defined bt the atoms H, N, CA and C of residues 3 to 5 of molecule 1.
   *
   * @class PeriodicTorsionProperty
   * @ingroup utils
   * @version Wed Jul 31 2002
   * @author gromos++ development team
   * @sa utils::Property
   */
  class PeriodicTorsionProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    PeriodicTorsionProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~PeriodicTorsionProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);
    /**
     * Calculate the torsional angle.
     */
    virtual Value const & calc();
    /**
     * calculate the nearest image distance
     */
    virtual Value nearestImageDistance(const Value & first, const Value & second) const;
  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
  };

  /**
   * Class CrossTorsionProperty
   * Purpose: Implements a cross torsion property.
   *
   * Description:
   * This class implements a cross torsion property. It is derived from the
   * Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim ct%<atomspec>%<step>%<lower_bound>%<upper_bound> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">\<step\></span>,
   *   <span style="color:darkred;font-family:monospace">\<lower_bound\></span> and
   *   <span style="color:darkred;font-family:monospace">\<upper_bound\></span>
   *   are additional parameters
   *
   * This calculates the torsion definied by the four atoms in the
   * <span style="color:darkred;font-family:monospace">\<atomspec\></span>
   * atom specifier. The periodic boundary conditions are taken into account.
   *
   * The additional parameters are only used by certain @ref programs "programs"
   * like @ref gca. Their meaning may change depdending on the program.
   *
   * For example:
   * - @verbatim t%1:1,2,3,4%1:3,4,5,6 @endverbatim means the cross torsion defined by the atoms 2,
   *   3, 4 and 5 and 3, 4, 5 and 6 of molecule 1.
   *
   * @class CrossTorsionProperty
   * @ingroup utils
   * @version 09.10.2009
   * @author N. Schmid
   * @sa utils::Property
   */

  class CrossTorsionProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    CrossTorsionProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~CrossTorsionProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Property parser mainly written for the HB Property
     */
    virtual void parse(AtomSpecifier const & atmspc);
    /**
     * Calculate the torsional angle.
     */
    virtual Value const & calc();
    /**
     * As most of the properties i can think of have something to do with
     * atoms and molecules, i define these vectors in the base class.
     * This is also done in order to be able to write one general
     * arguments parsing function in the base class.
     */
    AtomSpecifier & atoms2() { return d_atom2; }
    /**
     * and a const variant
     */
    AtomSpecifier const & atoms2()const { return d_atom2; }
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;
    /**
     * calculate the nearest image distance
     */
    virtual Value nearestImageDistance(const Value & first, const Value & second) const;
  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
    /**
     * atom specifier for the second dihedral
     */
    AtomSpecifier d_atom2;
  };

  /**
   * Class JValueProperty
   * Purpose: Implements a J-value property.
   *
   * Description:
   * This class implements a J-value property. It is derived from the 
   * Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim j%<atomspec>%<a>%<b>%<c>%<delta> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">\<atomspec\></span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">&lt;a&gt;</span>, 
   *   <span style="color:darkred;font-family:monospace">&lt;b&gt;</span> and
   *   <span style="color:darkred;font-family:monospace">&lt;c&gt;</span> are the
   *   parameters of the Karplus curve (defaults: 6.4, -1.4 and 1.9 Hz).
   * - <span style="color:darkred;font-family:monospace">&lt;delta&gt;</span> is an 
   *   angle @f$\delta@f$ in degree (default 0.0)
   *
   * This calculates the J-value for the torsion @f$\phi@f$ angle defined by 
   * the four atoms in the <span style="color:darkred;font-family:monospace">
   * &lt;atomspec&gt;</span> atom specifier by appling
   * @f[ J = a\cos(\phi + \delta)^2 + b\cos(\phi + \delta) + c\mathrm{,}@f]
   * where the angle @f$\delta@f$ is given by 
   * <span style="color:darkred;font-family:monospace">&lt;delta&gt;</span>. The 
   * parameters of the Karplus curve are given by 
   * <span style="color:darkred;font-family:monospace">&lt;a&gt;</span>, 
   * <span style="color:darkred;font-family:monospace">&lt;b&gt;</span> and
   * <span style="color:darkred;font-family:monospace">&lt;c&gt;</span>. Periodic
   * boundary conditions are taken into account.
   * 
   * For example:
   * - @verbatim j%1:1,3,4,5 @endverbatim means the J-value definied by the 
   *   torsion between the 1-3 and 4-5 bonds in molecule 1.
   * - @verbatim j%1:res(4:H,N,CA);va(1,1:res(4:CA,N,CB,C)) @endverbatim means 
   *   the J-value definied by the torsion between the H-N and CA-HA bonds of 
   *   residue 4 of molecule 1 where HA is a virtual atom defined by the atoms
   *   CA ,N, CB and C.
   *
   * Multiple J values are available via substitution.
   *
   * For example:
   * - @verbatim j%1:res((x):H,N,CA);va(1,1:res((x):CA,N,CB,C))|x=4,10 @endverbatim
   *   means the J values definied by the torsion between the H-N and CA-HA
   *   bonds of residue 4 and 10 of molecule 1  where HA is a virtual atom 
   *   defined by the atoms CA, N, CB and C.
   *
   * @class JValueProperty
   * @ingroup utils
   * @version Fri Dec 23 2005
   * @author gromos++ development team
   * @sa utils::Property
   */

  class JValueProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    JValueProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~JValueProperty();
    /**
     * Parse the arguments. Calls Property::parse.    * 
     */
       using Property::parse;
       virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculate the J-value of a torsional angle.
     */
    virtual Value const & calc();
  };    

  /**
   * Class VectorOrderProperty
   * Purpose: Implements a vector order property.
   *
   * Description:
   * Implements an order property (angle between two specified vectors)
   * It is derived from the Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim o%<vec1>%<vec2> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;vec1&gt;</span> and
   *   <span style="color:darkred;font-family:monospace">&lt;vec2&gt;</span> are
   *   @ref VectorSpecifier
   *
   * This calculates the order (angle) between two vectors defined by the
   * <span style="color:darkred;font-family:monospace">&lt;vec1&gt;</span> and 
   * <span style="color:darkred;font-family:monospace">&lt;vec2&gt;</span> vector
   * specifiers.
   *
   * For example:
   * - @verbatim o%cart(1,0,0)%atom(1:1,2) @endverbatim means the order between 
   *   the x axis and the vector definied by atoms 1 and 2 of the first molecule.
   * - @verbatim o%atom(1:res(2:H,N))%cart(0,0,1) @endverbatim means the order 
   *   between the N-H bond of residue 2 and the z axis.
   *
   * Multiple orders are available via substitution.
   *
   * For example:
   * - @verbatim o%atom(1:res((x):H,N))%cart(0,0,1)|x=4,8 @endverbatim means the
   *   orders between the N-H bonds of residues 4 and 8 of the first molecule
   *   and the z axis.
   *
   * @class VectorOrderProperty
   * @ingroup utils
   * @version Mar 22 2005
   * @author gromos++ development team
   * @sa utils::Property utils::VectorOrderParamProperty
   */
  class VectorOrderProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    VectorOrderProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~VectorOrderProperty();
    /**
     * Parse and check property specifier (given in arguments).
     * Calls Property::parse and checks the arguments.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculate the angle between the given atoms.
     */
    virtual Value const & calc();
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;

  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * vector 1
     */
    utils::VectorSpecifier d_vec1;
    /**
     * vector 1
     */
    utils::VectorSpecifier d_vec2;
    
  };

  /**
   * Class VectorOrderParamProperty
   * Purpose: Implements an order parameter property.
   *
   * Description:
   * Implements an order parameter property
   * (order parameter of the angle between two specified vectors)
   * It is derived from the Property class.
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim op%<vec1>%<vec2> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;vec1&gt;</span> and
   *   <span style="color:darkred;font-family:monospace">&lt;vec2&gt;</span> are 
   *   @ref VectorSpecifier
   *
   * This calculates the order parameter @f$o@f$ from the angle @f$\alpha@f$ 
   * between the two vectors specified by 
   * <span style="color:darkred;font-family:monospace">&lt;vec1&gt;</span> and
   * <span style="color:darkred;font-family:monospace">&lt;vec2&gt;</span> by
   * @f[o = \frac{1}{2}(3\cos(\alpha)^2 - 1)\mathrm{.}@f]
   *
   * For example:
   * - @verbatim op%cart(1,0,0)%atom(1:1,2) @endverbatim means the order
   *   parameter between the x axis and the vector definied by atoms 1 and 2 of 
   *   the first molecule.
   * - @verbatim op%atom(1:res(2:H,N))%cart(0,0,1) @endverbatim means the order
   *   parameter between the N-H bond of residue 2 and the z axis.
   *
   * Multiple orders are available via substitution.
   *
   * For example:
   * - @verbatim op%atom(1:res((x):H,N))%cart(0,0,1)|x=4,8 @endverbatim means
   *   the order parameters between the N-H bonds of residues 4 and 8 of the 
   *   first molecule and the z axis.
   *
   * @class VectorOrderParamProperty
   * @ingroup utils
   * @version Jan 16 2004
   * @author gromos++ development team
   * @sa utils::Property utils::VectorOrderProperty
   */

  class VectorOrderParamProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    VectorOrderParamProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~VectorOrderParamProperty();
    /**
     * Parse and check property specifier (given in arguments).
     * Calls Property::parse and checks the arguments.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculate the angle between the given atoms.
     */
    virtual Value const & calc();
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;

  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * vector 1
     */
    utils::VectorSpecifier d_vec1;
    /**
     * vector 1
     */
    utils::VectorSpecifier d_vec2;

  };

  /**
   * Class PseudoRotationProperty
   * Purpose: Implements a property that can calculate a pseudo rotation.
   *
   * Description:
   * This class implements a pseudo rotation. It is derived from the 
   * Property class. 
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim pr%<atomspec> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> is an
   *   @ref AtomSpecifier
   *
   * This calculates the pseudo rotation of five atoms given by the 
   * <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;></span>
   * atom specifier. The pseudo rotation (@f$\Delta / 2@f$) is defined according  
   * to:
   * Altona, C; Geise, HJ; Romers C; Tetrahedron 24 13-32 (1968) 
   *
   * For example:
   * - @verbatim pr%1:res(1:C1*,C2*,C3*,C4*,O4*) @endverbatim means the pseudo
   *   rotation of the atoms in the furanose ring of residue 1 of molecule 1.
   *
   * Multiple pseudo rotations are available via substitution.
   *
   * For example:
   * - @verbatim pr%1:res((x):C1*,C2*,C3*,C4*,O4*)|x=1-6 @endverbatim means the 
   *   pseudo rotations of the atoms in the furanose ring of residues 1 to 6 of
   *   molecule 1.
   *
   * @class PseudoRotationProperty
   * @ingroup utils
   * @version Fri Apr 23 2004
   * @author gromos++ development team
   * @sa utils::Property utils::PuckerAmplitudeProperty
   */

  class PseudoRotationProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    PseudoRotationProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~PseudoRotationProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculate the torsional angle.
     */
    virtual Value const & calc();

  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * Function to calculate a torsion for four atoms
     */
    double _calcDihedral(int const a, int const b, int const c, int const d);
    /**
     * A constant that is needed every time
     * Should be sin(36) + sin(72);
     */
    const double d_sin36sin72;
  };      

  /**
   * Class PuckerAmplitudeProperty
   * Purpose: Implements a property that can calculate a the amplitude of a 
   * pucker rotation.
   *
   * Description:
   * This class implements a pucker amplitude. It is derived from the 
   * PseudoRotationProperty class. 
   *
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim pa%<atomspec> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> is an
   *   @ref AtomSpecifier
   *
   * This calculates the pucker amplitude for the five atoms given by
   * &lt;atomspec&gt;. The amplitude is defined according  to:
   * Altona, C; Geise, HJ; Romers C; Tetrahedron 24 13-32 (1968) 
   * see also:
   * Altona, C; Sundaralingam, M; JACS 94 8205-8212 (1972)
   *
   * For example:
   * - @verbatim pa%1:res(1:C1*,C2*,C3*,C4*,O4*) @endverbatim means the pucker 
   *   amplitude of the atoms in the furanose ring of residue 1 of molecule 1.
   *
   * Multiple pseudo rotations are available via substitution.
   *
   * For example:
   * - @verbatim pa%1:res((x):C1*,C2*,C3*,C4*,O4*)|x=1-6 @endverbatim means the
   *   pucker aplitudes of the atoms in the furanose ring of residues 1 to 6 of 
   *   molecule 1.
   *
   * @class PuckerAmplitudeProperty
   * @ingroup utils
   * @version Fri Apr 23 2004
   * @author gromos++ development team
   * @sa utils::PseudoRotationProperty
   * @sa utils::Property
   */

  class PuckerAmplitudeProperty : public PseudoRotationProperty
  {
  public:
    /**
     * Constructor.
     */
    PuckerAmplitudeProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~PuckerAmplitudeProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual Value const & calc();
  };    
  
  /**
   * Class CPAmplitudeProperty
   * Purpose: Implements a property that can calculate the Cremer-Pople Amplitude Q.
   *
   * Description:
   * This class implements Cremer-Pople Amplitude Q for analysis of hexopyranoses. It is derived from the 
   * Property class. 
   *
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim cpq%<atomspec> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace"><atomspec></span> is an
   *   @ref AtomSpecifier
   *
   * This calculates the Cremer-Pople Amplitude Q for analysis of hexopyranoses. The 6 
   * ring atoms have to be specified in the correct order starting from the ring 
   * oxygen to be consistent with the definition of Cremer-Pople coordinates.
   *
   * For example:
   * - @verbatim cpq%1:res(1:O5*,C1*,C2*,C3*,C4*,C5*) @endverbatim means the Cremer
   *   Pople Q of a pyranose ring of residue 1 of molecule 1. 
   *
   * @class CPAmplitudeProperty
   * @ingroup utils
   * @version Thu, Aug 25th 2016
   * @author dahahn
   * @sa utils::CPAmplitudeProperty
   * @sa utils::Property
   */

  class CPAmplitudeProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    CPAmplitudeProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~CPAmplitudeProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    virtual Value const & calc();
  };

  /**
   * Class CPThetaProperty
   * Purpose: Implements a property that can calculate the Cremer-Pople Theta.
   *
   * Description:
   * This class implements Cremer-Pople Theta for analysis of hexopyranoses. It is derived from the 
   * Property class. 
   *
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim cpt%<atomspec> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace"><atomspec></span> is an
   *   @ref AtomSpecifier
   *
   * This calculates the Cremer-Pople Theta for analysis of hexopyranoses. The 6 
   * ring atoms have to be specified in the correct order starting from the ring 
   * oxygen to be consistent with the definition of Cremer-Pople coordinates.
   *
   * For example:
   * - @verbatim cpt%1:res(1:O5*,C1*,C2*,C3*,C4*,C5*) @endverbatim means the Cremer
   *   Pople theta of a pyranose ring of residue 1 of molecule 1. 
   *
   * @class CPThetaProperty
   * @ingroup utils
   * @version Tue, Aug 25th 2016
   * @author dahahn
   * @sa utils::CPThetaProperty
   * @sa utils::Property
   */

  class CPThetaProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    CPThetaProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~CPThetaProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    virtual Value const & calc();
  };

  /**
   * Class CPPhiProperty
   * Purpose: Implements a property that can calculate the Cremer-Pople Phi.
   *
   * Description:
   * This class implements Cremer-Pople Phi for analysis of hexopyranoses. It is derived from the 
   * Property class. 
   *
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim cpp%<atomspec> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace"><atomspec></span> is an
   *   @ref AtomSpecifier
   *
   * This calculates the Cremer-Pople Phi for analysis of hexopyranoses. The 6 
   * ring atoms have to be specified in the correct order starting from the ring 
   * oxygen to be consistent with the definition of Cremer-Pople coordinates.
   *
   * For example:
   * - @verbatim cpp%1:res(1:O5*,C1*,C2*,C3*,C4*,C5*) @endverbatim means the Cremer
   *   Pople phi of a pyranose ring of residue 1 of molecule 1. 
   *
   * @class CPPhiProperty
   * @ingroup utils
   * @version Tue, Jan 26th 2016
   * @author dahahn
   * @sa utils::CPPhiProperty
   * @sa utils::Property
   */

  class CPPhiProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    CPPhiProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~CPPhiProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    virtual void parse(std::vector<std::string> const & arguments, int x);
    virtual Value const & calc();
  };

  /**
   * Class ExpressionProperty
   * Purpose: Implements an expression property
   *
   * Description:
   * The expression property allows the evaluation of a multitude of expressions
   * over a trajectory.
   * The general form is:
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim expr%<f1>(<args1>) <op> <f2>(<args2>) @endverbatim
   * </b></span>
   * <br>
   * where:
   *  - <span style="color:darkred;font-family:monospace">&lt;op&gt;</span> is an
   *    operator:
   *    - arithmetic:  + - * /
   *    - logic: ! (not), ==, !=, <, >, <=, >=, &&, ||
   *    - <span style="color:darkred;font-family:monospace">&lt;f1&gt;</span> and 
   *      <span style="color:darkred;font-family:monospace">&lt;f2&gt;</span> are 
   *      functions. The functions take different kinds of arguments dependent on
   *      the function.
   *    - scalar as argument: sin, cos, tan, asin, acos, atan, exp, ln, abs, sqrt
   *    - vector as argument: abs, abs2 (squared abs)
   *    - two vectors as arguments: dot, cross, ni (nearest image)
   * - <span style="color:darkred;font-family:monospace">&lt;args1&gt;</span> and 
   *   <span style="color:darkred;font-family:monospace">&lt;args2&gt;</span> are 
   *   arguments for the functions and can be:
   *    - scalar
   *    - vector: defined by a @ref VectorSpecifier
   *
   * This calculates the expression and returns the scalar or the vector as 
   * result. Please note that spacing can influence the expression. Always, be
   * careful when using spaces and check whether the result was not affected.
   * 
   * For example:
   * - @verbatim expr%dot(atom(1:1),cart(0,0,1)) @endverbatim calculates the
   *   dot product between position of atom(1:1) and the vector (0,0,1); that
   *   is the z-component of the position of the first atom of the first 
   *   molecule.
   * - @verbatim expr%abs(atom(1:1) - ni(atom(2:1), atom(1:1))) @endverbatim
   *   calculates the distance between the first atom of the first and second 
   *   molecule. First, the nearest image of atom(2:1) according to atom(1:1) is
   *   calculated. Second, this vector is substracted from atom(1:1) and the
   *   absolute value is taken. This is the same as property 
   *   @verbatim d%1:1;2:1 @endverbatim
   * - @verbatim expr%abs(atom(1:1) - ni(atom(2:1), atom(1:1)))<1.0 @endverbatim
   *   returns 1 if the the discussed distance is below 1.0 nm and 0 if not.
   * - @verbatim expr%acos(dot(atom(1:1,2),cart(0,0,1)) / abs(atom(1:1,2)))*(180/3.1415) @endverbatim
   *   calculates the order of the vector defined by atoms 1:1,2 and the z axis.
   *   First, the cosine is calculated by the dot product of the vectors and 
   *   division by their lengths (the second vector has length 1). Them the 
   *   angle is calculated and converted to degrees. This corresponds to 
   *   @verbatim o%atom(1:1,2)%cart(0,0,1) @endverbatim
   *
   * Multiple expressions are available via substitution.
   *
   * For example:
   * - @verbatim expr%dot(atom(1:(x)),cart(0,0,1))|x=1,2,3 @endverbatim calculates
   *   the z component of atoms 1, 2 and 3 of molecule 1.
   *
   * @version Aug 25 2005
   * @author markus
   * @sa utils::Property
   */
  class ExpressionProperty : public Property
  {
  public:
    /**
     * Constructor.
     */
    ExpressionProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~ExpressionProperty();
    /**
     * Parse and check property specifier (given in arguments).
     * Calls Property::parse and checks the arguments.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculate the angle between the given atoms.
     */
    virtual Value const & calc();
    /**
     * Write a title string specifying the property.
     */
    virtual std::string toTitle()const;

  protected:
    /**
     * find the corresponding forcefield type of the property.
     * needs to be overwritten for the specific properties.
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);

    /**
     * calculates expressions
     */
    ExpressionParser<Value> d_parser;
    
    /**
     * stores the expression
     */
    std::vector<ExpressionParser<Value>::expr_struct> d_expr;
    
    /**
     * could store variables
     */
    std::map<std::string, Value> d_var;

  };

  /**
   * Class HBProperty
   * Purpose: Implements a HB property.
   *
   * Description:
   * This class implements a HB property. It is derived from the 
   * Property class.
   *
   * <b>Hydrogen Bonds</b>
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim hb%<atomspec>%<dist_upper>%<angle_lower> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span> is
   *   a distance (default: 0.25 nm)
   * - <span style="color:darkred;font-family:monospace">&lt;angle_lower&gt;</span> is
   *   an angle in degrees (default: 135 degree)
   *
   * The hydrogen bond is definied by three atoms given by the
   * <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> atom
   * specifier: the donor, the hydrogen and the acceptor. The second atom has to
   * be a hydrogen atom (check by mass).
   * 
   * The hydrogen bond is considered to be present (value of 1) if the 
   * hydrogen-acceptor distance is less than 
   * <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span> and
   * the angle  definied by the donor, hydrogen and acceptor atoms is larger than
   * <span style="color:darkred;font-family:monospace">&lt;angle_lower&gt;</span>.
   *
   * For exmaple:
   * - @verbatim hb%1:res(3:N,H);1:res(5:O) @endverbatim means the hyrdogen bond 
   *   between the H atom of residue 3 and the O atom of residue 5 of the first
   *   molecule.
   *
   * Multiple hydrogen bonds are available via substitution.
   *
   * For example:
   * - @verbatim hb%1:res((x):N,H);1:res((x+2):O)|x=3,4 @endverbatim means the 
   *   hydrogen bonds between H of residue 3 and 4 and O of the second next
   *   residue (5 and 6) of molecule 1.
   *
   * <b>Three Centered Hydrogen Bond</b>
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim hb%<atomspec>%<dist_upper>%<angle_lower>%<angle_sum>%<angle_plane> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span> is an
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span> is
   *   a distance (default: 0.27 nm)
   * - <span style="color:darkred;font-family:monospace">&lt;angle_lower&gt;</span>,
   *   <span style="color:darkred;font-family:monospace">&lt;angle_sum&gt;</span> and
   *   <span style="color:darkred;font-family:monospace">&lt;angle_plane&gt;</span>
   *   are an angles in degrees. (defaults: 90, 340, 15 degrees)
   *
   * A three centered hydrogen bond is defined for a donor atom D, hydrogen
   * atom H and two acceptor atoms A1 and A2 given by the atom specifier
   * <span style="color:darkred;font-family:monospace">&lt;atomspec&gt;</span>. 
   * The second atoms has to be a hydrogen atom (checked by mass).
   *
   * A three centered hydrogen bond is considered to be present (value of 1) if 
   * - the distances H-A1 and H-A2 are lower than 
   *   <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span>.
   * - the angles D-H-A1 and D-H-A2 are larger 
   *   <span style="color:darkred;font-family:monospace">&lt;angle_lower&gt;</span>. 
   * - the sum of the angles D-H-A1, D-H-A2 and A1-H-A2 is larger than 
   *   <span style="color:darkred;font-family:monospace">&lt;angle_sum&gt;</span>.
   * - the dihedral angle defined by the planes through the atoms D-A1-A2 and H-A1-A2 
   *   is smaller than <span style="color:darkred;font-family:monospace">
   *   &lt;angle_plane&gt;</span>
   *
   * For example:
   * - @verbatim hb%1:res(3:N,H);1:res(5:O);1:res(6:O) @endverbatim means the
   *   three centred hydrogen bond between atom H of residue 3 and atoms O of
   *   residues 5 and 6 of molecule 1.
   *
   * @class HBProperty
   * @ingroup utils
   * @version Mon Oct 31 2005
   * @author gromos++ development team
   * @sa utils::Property
   */

  class HBProperty : public Property
    {
  public:
    /**
     * Constructor.
     */
    HBProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~HBProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculates the HB.
     */
    virtual Value const & calc();
   /**
    * Returns the value in string representation.
    */
    virtual std::string toTitle()const;

  protected:
    /**
    * find the type of the property.
    * is case of HB: return 1 if 3center or 0 for no 3c 
    */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
    // member variables
    /** HB distance (2 if 3 center)
     */
    DistanceProperty d1_hb;
    DistanceProperty d2_hb;
    /** HB angles (3 if 3 center)
     */
    AngleProperty a1_hb;
    AngleProperty a2_hb;
    AngleProperty a3_hb;
    /** HB improper (1 if 3 center)
     */
    TorsionProperty i1_hb;
    /** test whether it is a HB (and 3cHB)
     */
    bool ishb;
    bool is3c;  
    };

  /**
   * Class StackingProperty
   * Purpose: Implements a stacking property.
   *
   * Description:
   * This class implements a stacking property. It is derived from the 
   * Property class.
   *
   * <b>Hydrogen Bonds</b>
   * <br>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim st%<atomspec1>%<atomspec1>%<dist_upper>%<angle_upper> @endverbatim
   * </b></span>
   * <br>
   * where:
   * - <span style="color:darkred;font-family:monospace">&lt;atomspec1&gt;</span> and
   *   <span style="color:darkred;font-family:monospace">&lt;atomspec2&gt;</span> are
   *   @ref AtomSpecifier
   * - <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span>
   *   is a distance (default: 0.5 nm)
   * - <span style="color:darkred;font-family:monospace">&lt;angle_upper&gt;</span> is
   *   an angle in degrees (default: 30 degree)
   *
   * A stacking interaction is defined for two ring systems defined by the
   * atom specifiers <span style="color:darkred;font-family:monospace">
   * &lt;atomspec1&gt;</span> and <span style="color:darkred;font-family:monospace">
   * &lt;atomspec2&gt;</span>. The first three atoms in the atom specifiers define the
   * plane through the ring. For the centre of geometry calculation all atoms in
   * the atom specifiers are taken into account.
   * 
   * Ring systems were considered to stack if the distance between the centres
   * of geometry of the rings is less than a given distance 
   * <span style="color:darkred;font-family:monospace">&lt;dist_upper&gt;</span> and the
   * angle between the planes through the two rings is less than the specified
   * angle <span style="color:darkred;font-family:monospace">&lt;angle_upper&gt;</span>.
   *
   * For exmaple:
   * - @verbatim st%1:1-3%1:4-6 @endverbatim means the stacking between the two
   *   planes definied by atoms 1-3 and 4-6 of the first molecule.
   * - @verbatim st%1:res(44:CG,CE1,CE2,CD1,CD2,CZ)%2:res(1:N1,C5,N3,C6,C2) @endverbatim
   *   means the stacking between this HISB ring of residue 44 of molecule 1 and
   *   and the pyrimidine ring of residue 2 of of molecule 2.
   *
   * Multiple stackings are available via substitution.
   *
   * For example:
   * - @verbatim st%1:1-3%2:res((x),N1,C5,N3,C6,C2)|x=1,2 @endverbatim means the
   *   stackings between the ring defined by 1:1-3 and the pyrimidine rings of
   *   residues 1 and 2 of molecule 2.
   *
   * @class StackingProperty
   * @ingroup utils
   * @version Son Oct 26 2008
   * @author ns
   * @sa utils::Property
   */

  class StackingProperty : public Property {
  public:
    /**
     * Constructor.
     */
    StackingProperty(gcore::System &sys, bound::Boundary * pbc);
    /**
     * Destructor.
     */
    virtual ~StackingProperty();
    /**
     * Parse the arguments. Calls Property::parse.
     */
    using Property::parse;
    virtual void parse(std::vector<std::string> const & arguments, int x);
    /**
     * Calculates the stacking.
     */
    virtual Value const & calc();
    /**
     * Returns the value in string representation.
     */
    virtual std::string toTitle()const;

    /**
     * accessor to the first plane
     */
    AtomSpecifier & atoms1() {
      return d_atoms1;
    }

    /**
     * accessor to the first plane
     */
    AtomSpecifier const & atoms1()const {
      return d_atoms1;
    }

    /**
     * accessor to the second plane
     */
    AtomSpecifier & atoms2() {
      return d_atoms2;
    }

    /**
     * accessor to the second plane
     */
    AtomSpecifier const & atoms2()const {
      return d_atoms2;
    }

  protected:
    /**
     * unsed for stacking
     */
    virtual int findTopologyType(gcore::MoleculeTopology const &mol_topo);
    // member variables

    /**
     * the sets of atoms for the planes
     */
    AtomSpecifier d_atoms1;
    AtomSpecifier d_atoms2;
  };
}

#endif

