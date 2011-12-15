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
  
  /**
   * @class BS_Dimension needed for BS_Vector to hold the value
   */
  class BS_Dimension {
  public:
    /**
     * The value of the coordinate
     */
    double value;
    /**
     * The periodicity of the variable. A periodicity of 0 means none.
     */
    double periodicity;
  };
    
  /**
   * @class Implementation for an n-dimensional vector holding the coordinates
   * of the subspace
   */
  class BS_Vector : public std::vector<BS_Dimension> {
  public:
    /**
     * return the length squared
     * @return squared length
     */
    double abs2();
    /**
     * this - subtrahend = result
     * @param[in]    subtrahend
     * @param[inout] result
     */
    void minus(const BS_Vector &subtrahend, BS_Vector &result);
    /**
     * scale the vector by scalar
     * @param[in] scalar
     */
    void scale(const double scalar);
    /**
     * Normalize the vector to length 1
     * @return the original length
     */
    double normalize();
    /**
     * Set every entry to zero.
     */
    void nullify();
    /**
     * Multiply the vector by scalar
     * @param[in] scalar
     * @return the scaled vector
     */
    BS_Vector operator*(const double scalar);
    /**
     * Add two vectors together
     * @param[in] summand
     * @return the sum
     */
    BS_Vector operator +(const BS_Vector &summand);
    /**
     * Add summand to a vector
     * @param summand
     */
    void operator +=(const BS_Vector &summand);
    /**
     * The dot product with other:
     *      dot(self, other)
     * @param other vector
     * @return the dot product
     */
    double dot(const BS_Vector &other);
    /**
     * Create a BS_Vector from two <double> vectors; one containing the values,
     * the other the periodicities.
     * @param values
     * @param periodicities
     */
    void create(std::vector<double> &values, std::vector<double> &periodicities);
    /**
     * Creates an output for the Vector.
     * @return the output
     */
    std::string str();
  };
  
  /**
   * @class BS_Coordinate
   * 
   * Reimplements the LE_Coordinate class, but with reduced coordinates,
   * so that different LE_Coordinates can be mixed.
   */
  class BS_Coordinate {
  public:
    BS_Coordinate(){
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
     * create a copy of the LE coordinate
     * @return a pointer to the copy
     */
    //virtual BS_Coordinate* clone() const = 0;
    /**
     * get type identifier
     * @return and integer describing the type
     */
    virtual int getType() const = 0;
    /**
     * Return the periodicity (0.0 if none)
     * @return the periodicity
     */
    virtual double getPeriodicity() = 0;
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
    virtual void getInternalCoordinates(std::vector<BS_Dimension> &coord) const = 0;
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
                             //m_periodicity = m_red_2pi * m_rad2degree;
                             m_periodicity = 0;
                             closeToZero = closeToPeriode = false;
    }
    virtual ~BS_Dihedral() {}
    //virtual BS_Dihedral* clone() const;
    virtual int getType() const { return m_type; /*dihedral*/ }
    virtual double getPeriodicity() {return m_periodicity;}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(std::vector<BS_Dimension> &coord) const;
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
    double m_periodicity;
    math::Vec fi, fj, fk, fl;
    // Used for undoing the refolding
    bool closeToZero, closeToPeriode;
  };

}
#endif	/* BS_LE_COORD_H */

