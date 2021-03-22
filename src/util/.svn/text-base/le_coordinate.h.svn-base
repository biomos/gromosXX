/**
 * @file le_coordinate.h
 * local elevation coordinates
 */

#ifndef INCLUDED_LE_COORDINATE_H
#define	INCLUDED_LE_COORDINATE_H

namespace configuration {
  class Configuration;
}

namespace util {
  /**
   * @class LE_Coordinate
   * @ingroup util
   * @short a LE generalized coordinate like a dihedral angle.
   *
   * This abstract class defines the interace for a LE coordinate. The interface is
   * implemented by different types of collective variables such
   * as dihedrals or distances
   *
   * @sa util::LE_Dihedral_Coordinate
   */
  class LE_Coordinate {
  public:
    LE_Coordinate() {
    }
    virtual ~LE_Coordinate() {
    }
    /**
     * get the umbrella ID
     * @return the umbrella ID of the LE Coordinate
     */
    int umbrella_id() const { return m_umbrella_id; }
    /**
     * set the umbrella ID
     * @param[in] id the umberlla ID to set
     */
    void umbrella_id(int id) { m_umbrella_id = id; }
    /**
     * create a copy of the LE coordinate
     * @return a pointer to the copy
     */
    virtual LE_Coordinate* clone() const = 0;
    /**
     * get type identifier
     * @return and integer describting the type
     */
    virtual int get_type() const = 0;
    /**
     * calculate the value and setup the pointers needed
     * @param[inout] conf the configuration for which the value is calculated.
     *             it will also be used to store the force
     */
    virtual void calculate(configuration::Configuration & conf) = 0;
    /**
     * get the value of the local elevation coordinate.
     * For periodic coordinates it will try to move to coordinate
     * to grid range
     * @param[in] grid_min the minimum value of the grid
     * @param[in] grid_max the maximum value of the grid
     * @return the value of the collective variable on the grid
     */
    virtual double get_value(double grid_min, double grid_max) const = 0;
    /**
     * get the deviation from the value, usually this is just the difference
     * @param[in] grid_value the value on the gird to which the deviation is
     *                   calculated
     * @return the deviation from the given value
     */
    virtual double get_deviation(const double & grid_value) const = 0;
    /**
     * apply the umbrella potential, calculate and add the force. The force
     * is added to the corresponding atoms in the configuration given in
     * the call.
     * The chain rule is applied to calculate the force
     * @f[ \mathbf{f}_i = - \frac{\partial{}V}{\partial\mathbf{r}_i} =
     *    - \frac{\partial{}V}{\partial{}Q}\frac{\partial{}Q}{\partial\mathbf{r}_i} @f]
     * @param[in] deriv the derivative of the potential energy by the collective
     *              variable @f$ \frac{\partial{}V}{\partial{}Q} @f$
     */
    virtual void apply(double deriv) = 0;
    /**
     * convert to string
     */
    virtual std::string str() const { return "LE_Coordinate"; }
  private:
    int m_umbrella_id;
  };

  /**
   * @class LE_Dihedral_Coordinate
   * @ingroup util
   * implementation of the LE dihedral angle coordinate. In addition to
   * the ID it is constructed from 4 atom indices defining the dihedral
   * angle
   * For detailed documentation see @ref util::LE_Coordinate
   *
   * @sa util::LE_Coordinate
   */
  class LE_Dihedral_Coordinate : public LE_Coordinate {
  public:
    /**
     * construct a new dihedral coordinate from atom indices
     * @param[in] id the umbrella ID
     * @param[in] i index of atom i
     * @param[in] j index of atom j
     * @param[in] k index of atom k
     * @param[in] l index of atom l
     */
    LE_Dihedral_Coordinate(int id, unsigned int i, unsigned int j,
            unsigned int k, unsigned int l) :
    i(i), j(j), k(k), l(l) { umbrella_id(id); }
    ~LE_Dihedral_Coordinate() {}
    virtual LE_Coordinate* clone() const;
    virtual int get_type() const { return 1; /*dihedral*/ }
    virtual void calculate(configuration::Configuration & conf);
    virtual double get_value(double grid_min, double grid_max) const;
    virtual double get_deviation(const double & grid_value) const;
    virtual void apply(double deriv);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    unsigned int i, j, k, l;
    configuration::Configuration * m_conf;
    double phi;
    math::Vec fi, fj, fk, fl;
  };
  /**
   * @class LE_Distance_Coordinate
   * @ingroup util
   * implementation of the LE distance coordinate. In addition to
   * the ID it is constructed from 2 atom indices defining the distance
   * For detailed documentation see @ref util::LE_Coordinate
   *
   * @sa util::LE_Coordinate
   */
  class LE_Distance_Coordinate : public LE_Coordinate {
  public:
    /**
     * construct a new distance coordinate from atom indices
     * @param[in] id the umbrella ID
     * @param[in] i index of atom i
     * @param[in] j index of atom j
     */
    LE_Distance_Coordinate(int id, unsigned int i, unsigned int j) :
    i(i), j(j) { umbrella_id(id); }
    ~LE_Distance_Coordinate() {}
    virtual LE_Coordinate* clone() const;
    virtual int get_type() const { return 2; /*distance*/ }
//    virtual int get_type() const { return util::umbrella::vt_distance; /*distance*/ }
    virtual void calculate(configuration::Configuration & conf);
    virtual double get_value(double grid_min, double grid_max) const;
    virtual double get_deviation(const double & grid_value) const;
    virtual void apply(double deriv);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    unsigned int i, j;
    configuration::Configuration * m_conf;
    double dist;
    math::Vec fi, fj;
  };
  /**
   * @class LE_DistanceField_Coordinate
   * @ingroup util
   * implementation of the LE distance field coordinate. In addition to
   * the ID it is constructed according to the input in the DISTANCEFIELD block
   * For detailed documentation see @ref util::LE_Coordinate
   *
   * @sa util::LE_Coordinate
   */
  class LE_DistanceField_Coordinate : public LE_Coordinate {
  public:
    /**
     * construct a new distance field coordinate
     * @param[in] id the umbrella ID
     */
    LE_DistanceField_Coordinate(int id, topology::Topology &topo, simulation::Simulation &sim);
    ~LE_DistanceField_Coordinate() {}
    virtual LE_Coordinate* clone() const;
    virtual int get_type() const { return 6; /*distance field*/ }
//    virtual int get_type() const { return util::umbrella::vt_distance; /*distance*/ }
    virtual void calculate(configuration::Configuration & conf);
    virtual double get_value(double grid_min, double grid_max) const;
    virtual double get_deviation(const double & grid_value) const;
    virtual void apply(double deriv);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    util::Virtual_Atom va_i, va_j;
    topology::Topology * m_topo;
    simulation::Simulation * m_sim;
    configuration::Configuration * m_conf;
    
    double dist;
    math::Vec fi, fj;
  };
}
#endif	/* INCLUDED_LE_COORDINATE_H */

