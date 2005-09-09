/**
 * @file gmath.h
 * mathematical definitions.
 */

#ifndef INCLUDED_MATH_H
#define INCLUDED_MATH_H

// #include <blitz/blitz.h>
// #include <blitz/array.h>
// #include <blitz/tinyvec-et.h>
// #include <blitz/tinymat.h>


namespace math
{

  // using namespace blitz;
  // BZ_USING_NAMESPACE(blitz)

  /**
   * 3 dimensional vector.
   */
  // typedef blitz::TinyVector<double, 3U> Vec;
  // typedef double Vec[3];
  class Vec
  {
  private:
    double d_v[3];
  public:
    Vec() {}
    explicit Vec(double d) { d_v[0] = d_v[1] = d_v[2] = d; }
    Vec(Vec const & v) { d_v[0] = v(0); d_v[1] = v(1); d_v[2] = v(2); }
    Vec(double d1, double d2, double d3) { d_v[0] = d1; d_v[1] = d2; d_v[2] = d3; }
    
    double operator()(int i)const { assert(i>=0 && i<3); return d_v[i]; }
    double & operator()(int i) { assert(i >= 0 && i < 3); return d_v[i]; }
    Vec & operator=(double d) { d_v[0] = d_v[1] = d_v[2] = d; return *this; }
    Vec & operator+=(Vec const &v) { d_v[0] += v(0); d_v[1] += v(1); d_v[2] += v(2); return *this; }
    Vec & operator-=(Vec const &v) { d_v[0] -= v(0); d_v[1] -= v(1); d_v[2] -= v(2); return *this; }
    Vec & operator*=(double d) { d_v[0] *= d; d_v[1] *= d; d_v[2] *= d; return *this; }
    Vec & operator/=(double d) { d_v[0] /= d; d_v[1] /= d; d_v[2] /= d; return *this; }
  };
  
  /**
   * Array of 3D vectors.
   */
  // typedef blitz::Array<Vec, 1>         VArray;
  // typedef std::vector<Vec> VArray;
  class VArray : public std::vector<Vec>
  {
  public:
#ifndef __SUNPRO_CC
    VArray() : std::vector<Vec>::vector() {}
    VArray(size_t s) : std::vector<Vec>::vector(s) {}
    VArray(size_t s, Vec const &v) : std::vector<Vec>::vector(s, v) {}
#else
    VArray() : vector() {}
    VArray(size_t s) : vector(s) {}
    VArray(size_t s, Vec const &v) : vector(s, v) {}
#endif
    
    VArray & operator=(double d)
    {
      for(std::vector<Vec>::iterator it=begin(), to=end(); it!=to; ++it)
	*it = d;
      return *this;
    }
    
    Vec const & operator()(int i)const 
    {
      assert(i < int(size()));
      return operator[](i);
    }
    Vec & operator()(int i)
    {
      assert(i < int(size()));
      return operator[](i);
    }

    VArray & operator+=(VArray const &v)
    {
      assert(size() == v.size());
      std::vector<Vec>::const_iterator it2 = v.begin();
      for(std::vector<Vec>::iterator it=begin(), to=end();
	  it!=to; ++it, ++it2)
	*it += *it2;
      return *this;
    }
  };
  
  /**
   * Array of scalars.
   */
  // typedef blitz::Array<double, 1>      SArray;
  // typedef std::vector<double> SArray;
  class SArray : public std::vector<double>
  {
  public:

#ifndef __SUNPRO_CC
    SArray() {}
    SArray(size_t s) : std::vector<double>::vector(s) {}
    SArray(size_t s, double d) : std::vector<double>::vector(s, d) {}
#else
    SArray() {}
    SArray(size_t s) : vector(s) {}
    SArray(size_t s, double d) : vector(s, d) {}
#endif
    double operator()(int i)const 
    {
      assert(i < int(size()));
      return operator[](i);
    }
    double & operator()(int i)
    {
      assert(i < int(size()));
      return operator[](i);
    }
  };
  
  /**
   * Matrix.
   */
  // typedef blitz::TinyMatrix<double, 3U, 3U> Matrix;
  // typedef double Matrix[3][3];
  class Matrix
  {
  private:
    double m[3][3];
  public:
    Matrix() {}
    explicit Matrix(double d) 
    {
      m[0][0] = m[0][1] = m[0][2] = 
	m[1][0] = m[1][1] = m[1][2] =
	m[2][0] = m[2][1] = m[2][2] = d;
    }
    
    Matrix & operator=(double d)
    {
      m[0][0] = m[0][1] = m[0][2] = 
	m[1][0] = m[1][1] = m[1][2] =
	m[2][0] = m[2][1] = m[2][2] = d;
      return *this;
    }
    
    Matrix & operator=(Matrix mat)
    {
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  m[i][j] = mat(i,j);
      return *this;
    }
    
    double operator()(int i, int j)const 
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return m[i][j];
    }
    double & operator()(int i, int j)
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return m[i][j];
    }
    Matrix & operator+=(Matrix const & mat)
    {
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  m[i][j] += mat(i,j);
      return *this;
    }
    Matrix & operator-=(Matrix const & mat)
    {
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  m[i][j] -= mat(i,j);
      return *this;
    }
  };

  /**
   * Box.
   */
  class Box
  {
  private:
    Vec d_b[3];
  public:
    Box() {}
    explicit Box(double d) 
    {
      d_b[0] = d;
      d_b[1] = d;
      d_b[2] = d;
    }
    Box(Vec const &v1, Vec const &v2, Vec const &v3)
    {
      d_b[0] = v1;
      d_b[1] = v2;
      d_b[2] = v3;
    }
    Vec const & operator()(int i)const { return d_b[i]; }
    Vec & operator()(int i) { return d_b[i]; }
    Box & operator*=(double d) { d_b[0] *= d; d_b[1] *= d; d_b[2] *= d; return *this; }

  };

  /**
   * @enum boundary_enum
   * boundary condition
   */
  enum boundary_enum{
    /**
     * vacuum.
     */
    vacuum = 0,
    /**
     * rectangular box
     */
    rectangular = 1,
    /**
     * triclinic box
     */
    triclinic = 2,
    /**
     * truncated octahedral box
     */
    truncoct = 3,
    /**
     * non-specialized version
     */
    any = 4
  };
  /**
   * @enum virial_enum
   * virial enum
   */
  enum virial_enum { 
    /** no virial */ 
    no_virial = 0, 
    /** molecular virial */
    molecular_virial = 1, 
    /** atomic virial */
    atomic_virial = 2 
  };
  /**
   * @enum pressure_scale_enum
   * pressure scaling
   */
  enum pressure_scale_enum{
    /** no pressure scaling */
    pcouple_off = 0,
    /** isotropic pressure scaling */
    pcouple_isotropic = 1,
    /** anisotropic pressure scaling */
    pcouple_anisotropic = 2,
    /** full anisotropic pressure scaling */
    pcouple_full_anisotropic = 3
  };

  /**
   * a small number.
   */
  const double epsilon = 0.000000000001;

  /**
   * Pi
   */
  const double Pi = 3.1415926535897932384626433;

  /**
   * Boltzmann constant.
   */
  extern double k_Boltzmann;

  /**
   * h bar.
   */
  extern double h_bar;

  /**
   * 1 / (4 Pi epsilon0).
   */
  extern double four_pi_eps_i;
  
#ifndef NDEBUG
  /**
   * module debug level.
   */
  extern int debug_level;
#endif


  inline bool operator==(Vec const &v1, Vec const &v2)
  {
    if (v1(0) != v2(0)) return false;
    if (v1(1) != v2(1)) return false;
    if (v1(2) != v2(2)) return false;
    return true;
  }

  /**
   * != operator
   */
  inline bool operator!=(Vec const &t1, Vec const &t2)
  {
    return !(t1 == t2);
  }

  inline double dot(Vec const &v1, Vec const &v2)
  {
    return v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2);
  }

  inline Vec cross(Vec const &v1, Vec const &v2)
  {
    return Vec(v1(1) * v2(2) - v1(2) * v2(1),
	       v1(2) * v2(0) - v1(0) * v2(2),
	       v1(0) * v2(1) - v1(1) * v2(0));
  }

  inline double abs2(Vec const &v)
  {
    return v(0) * v(0) + v(1) * v(1) + v(2) * v(2);
  }

  inline double abs(Vec const &v)
  {
    return sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2));
  }
  
  inline void dyade(Vec const &v1, Vec const &v2, Matrix &m)
  {
    for(int d1=0; d1 < 3; ++d1)
      for(int d2=0; d2 < 3; ++d2)
	m(d1,d2) = v1(d1) * v2(d2);
  }

  inline Vec operator+(Vec const &v1, Vec const &v2)
  {
    return Vec(v1(0) + v2(0),
	       v1(1) + v2(1),
	       v1(2) + v2(2));
  }

  inline Vec operator-(Vec const &v1)
  {
    return Vec(-v1(0),
	       -v1(1),
	       -v1(2));
  }
  
  inline Vec operator-(Vec const &v1, Vec const &v2)
  {
    return Vec(v1(0) - v2(0),
	       v1(1) - v2(1),
	       v1(2) - v2(2));
  }

  inline Vec operator*(Vec const &v1, double d)
  {
    return Vec(v1(0) * d,
	       v1(1) * d,
	       v1(2) * d);
  }

  inline Vec operator*(double d, Vec const &v1)
  {
    return Vec(v1(0) * d,
	       v1(1) * d,
	       v1(2) * d);
  }

  inline Vec operator/(Vec const &v1, double d)
  {
    return Vec(v1(0) / d,
	       v1(1) / d,
	       v1(2) / d);
  }

  inline std::string v2s(Vec const & v)
  {
    std::stringstream s;
    s << "[" 
      << std::setprecision(9) << std::setw(20) << v(0)
      << std::setprecision(9) << std::setw(20) << v(1)
      << std::setprecision(9) << std::setw(20) << v(2)
      << "]";
    return s.str();
  }

  inline double sum(SArray const & a)
  {
    double d = 0.0;
    for(unsigned int i=0; i<a.size(); ++i)
      d += a[i];
    return d;
  }

  inline Vec product(Matrix const &m, Vec const &v)
  {
    return Vec(
	       m(0,0) * v(0) + m(1,0) * v(1) + m(2,0) * v(2),
	       m(0,1) * v(0) + m(1,1) * v(1) + m(2,1) * v(2),
	       m(0,2) * v(0) + m(1,2) * v(1) + m(2,2) * v(2)
	       );
  }

  inline Matrix product(Matrix const &m1, Matrix const &m2)
  {
    Matrix m(0.0);
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	for(int k=0; k<3; ++k)
	  m(i,j) += m1(i,k) * m2(k,j);
    return m;
  }

  inline Box product(Matrix const &m1, Box const &m2)
  {
    Box m(0.0);
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	for(int k=0; k<3; ++k)
	  m(i)(j) += m1(i,k) * m2(k)(j);
    return m;
  }

  inline std::string m2s(Matrix const & m)
  {
    std::stringstream s;
    s << "[[" 
      << std::setprecision(9) << std::setw(20) << m(0,0)
      << std::setprecision(9) << std::setw(20) << m(0,1)
      << std::setprecision(9) << std::setw(20) << m(0,2)
      << "]\n"
      << "\t[" 
      << std::setprecision(9) << std::setw(20) << m(1,0)
      << std::setprecision(9) << std::setw(20) << m(1,1)
      << std::setprecision(9) << std::setw(20) << m(1,2)
      << "]\n"
      << "\t[" 
      << std::setprecision(9) << std::setw(20) << m(2,0)
      << std::setprecision(9) << std::setw(20) << m(2,1)
      << std::setprecision(9) << std::setw(20) << m(2,2)
      << "]]";
    
    return s.str();
  }


} // math

#endif

