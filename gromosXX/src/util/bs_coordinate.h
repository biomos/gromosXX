/**
 * @file:   bs_coordinate.h
 * Reimplements the LE_Coordinate classes, but with reduced coordinates
 */

#ifndef BS_COORDINATE_H
#define	BS_COORDINATE_H

namespace configuration {
  class Configuration;
}

namespace util{
  
  class BS_Vector;
  class BS_Dimension;
  /**
   * @class BS_Coordinate
   * 
   * Reimplements the LE_Coordinate class, but with reduced coordinates,
   * so that different LE_Coordinates can be mixed.
   */
  class BS_Coordinate {
  public:
    BS_Coordinate(){
      m_type = unknown;
    }
    virtual ~BS_Coordinate(){}
    /**
     * get the umbrella ID
     * @return the umbrella ID of the LE Coordinate
     */
    int umbrella_id() const { return m_umbrella_id; }
    /**
     * set the umbrella ID
     * @param[in] id the umbrella ID to set
     */
    void umbrella_id(int id) { m_umbrella_id = id; }
    /**
     * set the reduction factor for this coordinate
     * @param[in] rf The reduction factor
     */
    void red_fac(double rf) {m_red_fac = rf;}
    /**
     * Obtain the reduction factor
     * @return The reduction factor
     */
    double red_fac() {return m_red_fac;}
    /**
     * get type identifier
     * @return and integer describing the type
     */
    int getType() {return m_type;}
    /**
     * @return returns the dimensionality of the coordinate
     */
    int getDimension() {return m_dimension;};
    /**
     * calculate the internal Coordinate and the derivatives
     * @param[inout] conf the configuration for which the value is calculated.
     *             it will also be used to store the force
     */
    virtual void calculateInternalCoord(configuration::Configuration &conf) = 0;
    /**
     * Obtain the Internal Coordinates
     */
    virtual void getInternalCoordinates(BS_Vector &coord) const = 0;
    //virtual void getInternalCoordinates(std::vector<BS_Dimension> &coord) const = 0;
    /**
     * Add the forces expressed in internal Coordinates 
     * to the particles in conf.
     * @param [inout] conf the configuration to which the forces should be added
     * @param [inout] icForces the forces expressed in the particular internal
     *                         coordinates
     */
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives) = 0;
    /**
     * The enumeration defining, what kind of coordinate, we have
     */
    enum Coord_type {
        unknown = 0,
        dihedral,
        distance,
        lambda,
    };
    Coord_type m_type;
    /**
     * convert to string
     */
    virtual std::string str() const { return "BS_LE_Coordinate"; }

  protected:
    int m_umbrella_id;
    /**
     * The reduction factor for the coordinates
     */
    double m_red_fac;
    /**
     * The dimension of the internal Coordinates
     */
    int m_dimension;
    /**
     * The periodicity of the coordinate
     */
    //double m_periodicity;
 };
  
  /**
   * The BS_LE_Coordinate for dihedral angle. The internal coordinates are 
   * handled in degrees.
   * @class BS_LE_Dihedral
   *
   */
  class BS_Dihedral : public BS_Coordinate {
  public:
    /**
     * construct a new dihedral coordinate from atom indices
     * @param[in] id the umbrella ID
     * @param[in] i index of atom i
     * @param[in] j index of atom j
     * @param[in] k index of atom k
     * @param[in] l index of atom l
     * @param[in] red_fac reduction factor of this coordinate
     */
    BS_Dihedral(int id, unsigned int i, unsigned int j,
            unsigned int k, unsigned int l, double red_fac) :
    i(i), j(j), k(k), l(l) { umbrella_id(id); 
                             m_type = dihedral;
                             m_dimension = 1;
                             m_red_fac = red_fac;
                             m_red_pi = math::Pi / m_red_fac;
                             m_red_2pi = 2 * m_red_pi;
                             m_rad2degree = 180 / math::Pi;
                             m_counter = 0;
                             old_phi = math::Pi;
    }
    virtual ~BS_Dihedral() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    unsigned int i, j, k, l;
    double phi;
    double red_phi;
    double m_red_pi;
    double m_red_2pi;
    double m_rad2degree;
    math::Vec fi, fj, fk, fl;
    // Used for undoing the refolding
    double old_phi;
    int m_counter;
  };

  /**
   * @class BS_Distance
   * A internal Coordinate defined as distance between two atoms.
   */
  class BS_Distance : public BS_Coordinate {
  public:
    /**
     * Construct a new distance coordinate between atoma i and j.
     * @param id        The ID
     * @param i         The first atom
     * @param j         The second atom
     * @param red_fac   The reduction factor sigma
     */
    BS_Distance(int id, unsigned int i, unsigned j, double red_fac) :
    i(i), j(j) {
      umbrella_id(id);
      m_type = distance;
      m_red_fac = red_fac;
      m_dimension = 1;
    }
    virtual ~BS_Distance() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    unsigned int i, j;
    double m_distance;
    double m_reducedDistance;
    math::Vec fi, fj;
    
  };
}
#endif	/* BS_LE_COORD_H */

