/**
 * @file:   bs_coordinate.h
 * Reimplements the LE_Coordinate classes, but with reduced coordinates
 */

#ifndef BS_COORDINATE_H
#define	BS_COORDINATE_H

namespace configuration {
  class Configuration;
}

#include "bs_vector.h"

namespace util{
  
  /**
   * @class BS_Coordinate
   * @ingroup util
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
     * get the coordinate ID
     * @return the coordinate ID of the LE Coordinate
     */
    int cid() const { return m_id; }
    /**
     * set the coordinate ID
     * @param[in] id the coordinate ID to set
     */
    void cid(int id) { m_id = id; }
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
    unsigned int getDimension() {return m_dimension;};
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
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to) = 0;
    /**
     * Add the forces expressed in internal Coordinates.
     * The derivatives are dB / dQ.
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
        dihedral = 1,
        distance = 2,
        cartesian = 3,
        dihedralSum = 4,
        lambda = 5,
    };
    Coord_type m_type;
    /**
     * convert to string
     */
    virtual std::string str() const { return "BS_LE_Coordinate"; }
    /**
     * String, which is shown in the output file upon initialization
     */
    std::string init_str();

  protected:
    int m_id;
    /**
     * The reduction factor for the coordinates
     */
    double m_red_fac;
    /**
     * The dimension of the internal Coordinates
     */
    unsigned int m_dimension;
 };
  
  /**
   * @class BS_LE_Dihedral
   * 
   * The BS_LE_Coordinate for dihedral angle. The internal coordinates are 
   * handled in degrees.
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
            unsigned int k, unsigned int l, double red_fac);
    virtual ~BS_Dihedral() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to);
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
    double getPhi() const {return red_phi * m_rad2degree;}
    void incrCounter() {m_counter++;};
    void decrCounter() {m_counter--;};
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
    BS_Distance(int id, unsigned int i, unsigned int j, double red_fac) :
    i(i), j(j) {
      cid(id);
      m_type = distance;
      m_red_fac = red_fac;
      m_dimension = 1;
    }
    virtual ~BS_Distance() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to){}
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
  
  /**
   * @class BS_DihedralSum
   * 
   * The sum of two dihedral angles
   */
  class BS_DihedralSum : public BS_Coordinate {
  public:
    /**
     * The Sum of two dihedral angles
     * @param id    The ID
     * @param i     atom i of the first dihedral angle
     * @param j     atom j of the first dihedral angle
     * @param k     atom k of the first dihedral angle
     * @param l     atom l of the first dihedral angle
     * @param ii    atom i of the second dihedral angle
     * @param jj    atom j of the second dihedral angle
     * @param kk    atom k of the second dihedral angle
     * @param ll    atom l of the second dihedral angle
     * @param red_fac   The reduction factor
     */
    BS_DihedralSum(int id, unsigned int i, unsigned int j, unsigned int k,
            unsigned int l, unsigned int ii, unsigned int jj, unsigned int kk,
            unsigned int ll, double red_fac);
    virtual ~BS_DihedralSum() {};
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to);
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
    
  protected:
    unsigned int i, j, k, l, ii, jj, kk, ll;
    BS_Dihedral m_phi, m_psi;
    double m_sum;
    double m_old_sum;
    bool m_first;
  };
  
  class BS_Cartesian : public BS_Coordinate {
  public:
    BS_Cartesian(int id, std::vector<unsigned int> &atoms, 
                 bool allAtoms, double red_fac) {
      m_id = id;
      m_dimension = 3 * atoms.size();
      m_atoms = atoms;
      m_allAtoms = allAtoms;
      m_red_fac = red_fac;
      m_type = cartesian;
    }
    virtual ~BS_Cartesian() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to){}
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    std::vector<unsigned int> m_atoms;
    BS_Vector m_coordinates;
    bool m_allAtoms;
  };
  
  class BS_Lambda : public BS_Coordinate {
  public:
    BS_Lambda(int id, double red_fac) {
      m_id = id;
      m_dimension = 1;
      m_red_fac = red_fac;
      m_type = lambda;
    }
    virtual ~BS_Lambda() {}
    virtual void calculateInternalCoord(configuration::Configuration & conf);
    virtual void getInternalCoordinates(BS_Vector &coord) const;
    virtual void setOldPos(BS_Vector::iterator it, BS_Vector::iterator to){}
    virtual void addForces(configuration::Configuration &conf, 
                      BS_Vector &derivatives);
    virtual std::string str() const;
  protected:
    template<math::boundary_enum B>
    void _calculate(configuration::Configuration & conf);
    double m_lambda;
  };
}
#endif	/* BS_LE_COORD_H */

