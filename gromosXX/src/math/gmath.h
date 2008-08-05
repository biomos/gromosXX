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
  /**
   * 3 dimensional vector.
   */
  template<typename numeric_type>
  class GenericVec
  {
  private:
    numeric_type d_v[3];
  public:
    GenericVec() {}
    explicit GenericVec(numeric_type d) { d_v[0] = d_v[1] = d_v[2] = d; }
    template<typename numeric_type_b>
    GenericVec(GenericVec<numeric_type_b> const & v) { d_v[0] = v(0); d_v[1] = v(1); d_v[2] = v(2); }
    GenericVec(numeric_type d1, numeric_type d2, numeric_type d3) { d_v[0] = d1; d_v[1] = d2; d_v[2] = d3; }

    numeric_type operator()(int i)const { assert(i>=0 && i<3); return d_v[i]; }
    numeric_type & operator()(int i) { assert(i >= 0 && i < 3); return d_v[i]; }
    template<typename numeric_type_b>
    GenericVec<numeric_type> & operator=(numeric_type_b d) { d_v[0] = d_v[1] = d_v[2] = d; return *this; }
    template<typename numeric_type_b>
    GenericVec<numeric_type> & operator+=(GenericVec<numeric_type_b> const &v) { d_v[0] += v(0); d_v[1] += v(1); d_v[2] += v(2); return *this; }
    template<typename numeric_type_b>
    GenericVec<numeric_type> & operator-=(GenericVec<numeric_type_b> const &v) { d_v[0] -= v(0); d_v[1] -= v(1); d_v[2] -= v(2); return *this; } 
    template<typename numeric_type_b>
    GenericVec<numeric_type> & operator*=(numeric_type_b d) { d_v[0] *= d; d_v[1] *= d; d_v[2] *= d; return *this; }
    template<typename numeric_type_b>
    GenericVec<numeric_type> & operator/=(numeric_type_b d) { d_v[0] /= d; d_v[1] /= d; d_v[2] /= d; return *this; }
    };
    
    /**
     * double vector
     */
    typedef GenericVec<double> Vec;
    /**
     * long double vector
     */
    typedef GenericVec<long double> Vecl;
  
  /**
   * Array of 3D vectors.
   */
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
      assert(i >= 0 && i < int(size()));
      return operator[](i);
    }
    Vec & operator()(int i)
    {
      assert(i >= 0 && i < int(size()));
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
      assert(i >= 0 && i < int(size()));
      return operator[](i);
    }
    double & operator()(int i)
    {
      assert(i >= 0 && i < int(size()));
      return operator[](i);
    }
  };
  
  // the box should be known to the Matrix
  class Box;
  
  /**
   * Matrix.
   */
  template<typename numeric_type>
  class GenericMatrix
  {
  private:
    numeric_type m[3][3];
  public:
    GenericMatrix() {}
    template<typename numeric_type_b>
    explicit GenericMatrix(numeric_type_b d) 
    {
      m[0][0] = m[0][1] = m[0][2] = 
	m[1][0] = m[1][1] = m[1][2] =
	m[2][0] = m[2][1] = m[2][2] = d;
    }

    template<typename numeric_type_b>
    GenericMatrix(const GenericVec<numeric_type_b> &u,
           const GenericVec<numeric_type_b> &v,
           const GenericVec<numeric_type_b> &w,
           bool column_wise = false) {
      if (column_wise) {
        for (int i = 0; i < 3; ++i) {
          m[i][0] = u(i);
          m[i][1] = v(i);
          m[i][2] = w(i);
        }
      } else {
        for (int i = 0; i < 3; ++i) {
          m[0][i] = u(i);
          m[1][i] = v(i);
          m[2][i] = w(i);
        }
      }
    }
    
    GenericMatrix(const Box & box);

    template<typename numeric_type_b>
    inline GenericMatrix(
            const numeric_type_b &d1, const numeric_type_b &d2, const numeric_type_b &d3,
            const numeric_type_b &d4, const numeric_type_b &d5, const numeric_type_b &d6,
            const numeric_type_b &d7, const numeric_type_b &d8, const numeric_type_b &d9) {
      m[0][0] = d1;
      m[0][1] = d2;
      m[0][2] = d3;
      m[1][0] = d4;
      m[1][1] = d5;
      m[1][2] = d6;
      m[2][0] = d7;
      m[2][1] = d8;
      m[2][2] = d9;
    }
   
    template<typename numeric_type_b>
    GenericMatrix & operator=(numeric_type_b d)
    {
      m[0][0] = m[0][1] = m[0][2] = 
	m[1][0] = m[1][1] = m[1][2] =
	m[2][0] = m[2][1] = m[2][2] = d;
      return *this;
    }
    
    template<typename numeric_type_b>
    GenericMatrix & operator=(GenericMatrix<numeric_type_b> mat)
    {
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  m[i][j] = mat(i,j);
      return *this;
    }
    
    numeric_type operator()(int i, int j)const 
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return m[i][j];
    }
    numeric_type & operator()(int i, int j)
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return m[i][j];
    }

    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator+=(GenericMatrix<numeric_type_b> const & mat) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          m[i][j] += mat(i, j);
      return *this;
    }

    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator-=(GenericMatrix<numeric_type_b> const & mat) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          m[i][j] -= mat(i, j);
      return *this;
    }
    
    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator*=(const numeric_type_b & d) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          m[i][j] *= d;
      return *this;
    }
    
    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator/=(const numeric_type_b & d) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          m[i][j] /= d;
      return *this;
    }
  };
  
  /**
   * a double matrix
   */
  typedef GenericMatrix<double> Matrix;
  /**
   * a long double matrix
   */
  typedef GenericMatrix<long double> Matrixl;

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
    double operator()(int i, int j)const 
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return d_b[j](i);
    }
    Vec & operator()(int i) { return d_b[i]; }
    Box & operator*=(double d) { d_b[0] *= d; d_b[1] *= d; d_b[2] *= d; return *this; }
    Box & operator/=(double d) { d_b[0] /= d; d_b[1] /= d; d_b[2] /= d; return *this; }
    Box & operator=(Box const & box) 
    {
      d_b[0] = box.d_b[0];
      d_b[1] = box.d_b[1];
      d_b[2] = box.d_b[2];
      return *this;
    }
    inline Box & operator+=(Box const & box) {
      d_b[0] += box.d_b[0];
      d_b[1] += box.d_b[1];
      d_b[2] += box.d_b[2];
      return *this;
    }
    inline Box & operator-=(Box const & box) {
      d_b[0] -= box.d_b[0];
      d_b[1] -= box.d_b[1];
      d_b[2] -= box.d_b[2];
      return *this;
    }
  };

  template<typename numeric_type>
  GenericMatrix<numeric_type>::GenericMatrix(const Box & box) {
    for (unsigned int i = 0; i < 3; ++i) {
      m[i][0] = box(0)(i);
      m[i][1] = box(1)(i);
      m[i][2] = box(2)(i);
    }
  }

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
  /**
   * @f$ \epsilon^{-1} @f$
   */
  extern double eps0_i;
  
#ifndef NDEBUG
  /**
   * module debug level.
   */
  extern int debug_level;
#endif

  template<typename numeric_type_a, typename numeric_type_b>
  inline bool operator==(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    if (v1(0) != v2(0)) return false;
    if (v1(1) != v2(1)) return false;
    if (v1(2) != v2(2)) return false;
    return true;
  }

  /**
   * != operator
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline bool operator!=(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    return !(v1 == v2);
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline double dot(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    return v1(0) * v2(0) + v1(1) * v2(1) + v1(2) * v2(2);
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> cross(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    return GenericVec<numeric_type_a>(v1(1) * v2(2) - v1(2) * v2(1),
	       v1(2) * v2(0) - v1(0) * v2(2),
	       v1(0) * v2(1) - v1(1) * v2(0));
  }

  template<typename numeric_type>
  inline double abs2(GenericVec<numeric_type> const &v)
  {
    return v(0) * v(0) + v(1) * v(1) + v(2) * v(2);
  }

  template<typename numeric_type>
  inline double abs(GenericVec<numeric_type> const &v)
  {
    return sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2));
  }
  
  inline void dyade(Vec const &v1, Vec const &v2, Matrix &m)
  {
    for(int d1=0; d1 < 3; ++d1)
      for(int d2=0; d2 < 3; ++d2)
	m(d1,d2) = v1(d1) * v2(d2);
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> operator+(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    return Vec(v1(0) + v2(0),
	       v1(1) + v2(1),
	       v1(2) + v2(2));
  }

  template<typename numeric_type>
  inline GenericVec<numeric_type> operator-(GenericVec<numeric_type> const &v1)
  {
    return Vec(-v1(0),
	       -v1(1),
	       -v1(2));
  }
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> operator-(GenericVec<numeric_type_a> const &v1, GenericVec<numeric_type_b> const &v2)
  {
    return Vec(v1(0) - v2(0),
	       v1(1) - v2(1),
	       v1(2) - v2(2));
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> operator*(GenericVec<numeric_type_a> const &v1, numeric_type_b d)
  {
    return Vec(v1(0) * d,
	       v1(1) * d,
	       v1(2) * d);
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> operator*(numeric_type_b d, GenericVec<numeric_type_a> const &v1)
  {
    return Vec(v1(0) * d,
	       v1(1) * d,
	       v1(2) * d);
  }

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_a> operator/(GenericVec<numeric_type_a> const &v1, numeric_type_b d)
  {
    return Vec(v1(0) / d,
	       v1(1) / d,
	       v1(2) / d);
  }

  inline std::string v2s(GenericVec<double> const & v)
  {
    std::stringstream s;
    s << "[" 
      << std::setprecision(9) << std::setw(20) << v(0)
      << std::setprecision(9) << std::setw(20) << v(1)
      << std::setprecision(9) << std::setw(20) << v(2)
      << "]";
    return s.str();
  }
  
  inline std::string v2s(GenericVec<int> const & v)
  {
    std::stringstream s;
    s << "[" 
      << std::setw(8) << v(0)
      << std::setw(8) << v(1)
      << std::setw(8) << v(2)
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

  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericVec<numeric_type_b> product(GenericMatrix<numeric_type_b> const &m, GenericVec<numeric_type_a> const &v)
  {
    return GenericVec<numeric_type_b>(
	       m(0,0) * v(0) + m(1,0) * v(1) + m(2,0) * v(2),
	       m(0,1) * v(0) + m(1,1) * v(1) + m(2,1) * v(2),
	       m(0,2) * v(0) + m(1,2) * v(1) + m(2,2) * v(2)
	       );
  }
  
  inline Vec product(Box const &m, Vec const &v)
  {
    return Vec(
	       m(0,0) * v(0) + m(1,0) * v(1) + m(2,0) * v(2),
	       m(0,1) * v(0) + m(1,1) * v(1) + m(2,1) * v(2),
	       m(0,2) * v(0) + m(1,2) * v(1) + m(2,2) * v(2)
	       );
  }
  
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_a> product(
          GenericMatrix<numeric_type_a> const &m1,
          GenericMatrix<numeric_type_b> const &m2)
  {  
    GenericMatrix<numeric_type_a> m(numeric_type_a(0.0));
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          m(i, j) += m1(i,k) * m2(k,j);
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
  
  inline Matrix product(Box const &m1, Matrix const &m2)
  {
    Matrix m(0.0);   
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	for(int k=0; k<3; ++k)
	  m(i,j) += m1(i)(k) * m2(k,j);
    return m;
  }
  
  inline std::string m2s(GenericMatrix<double> const & m)
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
  
  inline int sign (double signum ){
    if(signum<0)
        return -1;
    else
        return 1;
  }

  inline long double costest (double acos_param){
      if(fabs(acos_param)>1)
          if(fabs(acos_param)>1+epsilon)
              return sign(acos_param)*1.0; 
          else 
              return acos_param;
      else return acos_param;
  }
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_a> operator*(GenericMatrix<numeric_type_a> const &ma, numeric_type_b d)
  {
    return GenericMatrix<numeric_type_a>(ma(0,0) * d, ma(0,1) * d, ma(0,2) * d,
            ma(1,0) * d, ma(1,1) * d, ma(1,2) * d,
            ma(2,0) * d, ma(2,1) * d, ma(2,2) * d);
  }
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_a> operator+(GenericMatrix<numeric_type_a> const &ma, GenericMatrix<numeric_type_b> const &mb)
  {
    GenericMatrix<numeric_type_a> m;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        m(i,j) = ma(i,j) + mb(i,j);
    return m;
  }
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_a> operator-(GenericMatrix<numeric_type_a> const &ma, GenericMatrix<numeric_type_b> const &mb)
  {
    GenericMatrix<numeric_type_a> m;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        m(i,j) = ma(i,j) - mb(i,j);
    return m;
  }
  
  /**
   * square a matric
   */
  template<typename numeric_type>
  inline GenericMatrix<numeric_type> square(GenericMatrix<numeric_type> const &ma)
  {
    GenericMatrix<numeric_type> m;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        m(i,j) = ma(i,j) * ma(i,j);
    return m;
  } 
  
  /**
   * matrix to the power of x
   */
  template<typename numeric_type>
  inline GenericMatrix<numeric_type> square(GenericMatrix<numeric_type> const &ma, const double & x)
  {
    GenericMatrix<numeric_type> m;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        m(i,j) = std::pow(ma(i,j), x);
    return m;
  } 
  /**
   * determinat of A
   */
  template<typename numeric_type>
  inline numeric_type det(GenericMatrix<numeric_type> const &ma)
  {
    
    return ma(0,1)*ma(1,2)*ma(2,0)
          -ma(0,2)*ma(1,1)*ma(2,0)
          +ma(0,2)*ma(1,0)*ma(2,1)
          -ma(0,0)*ma(1,2)*ma(2,1)
          -ma(0,1)*ma(1,0)*ma(2,2)
          +ma(0,0)*ma(1,1)*ma(2,2);
  }
  
  /**
   * inverse of Matrix X
   */
  template<typename numeric_type>
  inline GenericMatrix<numeric_type> inverse(GenericMatrix<numeric_type> const &ma){
    GenericMatrix<numeric_type> inv;
    const double det_ma = det(ma);
    inv(0, 0) = (ma(1, 1) * ma(2, 2) - ma(1, 2) * ma(2, 1)) / det_ma;
    inv(1, 0) = (ma(1, 2) * ma(2, 0) - ma(1, 0) * ma(2, 2)) / det_ma;
    inv(2, 0) = (ma(1, 0) * ma(2, 1) - ma(1, 1) * ma(2, 0)) / det_ma;
    inv(0, 1) = (ma(0, 2) * ma(2, 1) - ma(0, 1) * ma(2, 2)) / det_ma;
    inv(1, 1) = (ma(0, 0) * ma(2, 2) - ma(0, 2) * ma(2, 0)) / det_ma;
    inv(2, 1) = (ma(0, 1) * ma(2, 0) - ma(0, 0) * ma(2, 1)) / det_ma;
    inv(0, 2) = (ma(0, 1) * ma(1, 2) - ma(0, 2) * ma(1, 1)) / det_ma;
    inv(1, 2) = (ma(0, 2) * ma(1, 0) - ma(0, 0) * ma(1, 2)) / det_ma;
    inv(2, 2) = (ma(0, 0) * ma(1, 1) - ma(0, 1) * ma(1, 0)) / det_ma;
    return inv;
  }
  
  /**
   * transpose of Matrix X
   */
  template<typename numeric_type>
  inline GenericMatrix<numeric_type> transpose(GenericMatrix<numeric_type> const &ma) {
    GenericMatrix<numeric_type> trans;
    for(unsigned int i = 0; i < 3; ++i) {
      for(unsigned int j = 0; j < 3; ++j) {
        trans(j,i) = ma(i, j);
      }
    }
    return trans;
  }
 
  /**
   * trace of Matrix X
   */
  template<typename numeric_type>
  inline numeric_type trace(GenericMatrix<numeric_type> const &ma) {
    return ma(0,0) * ma(1,1) * ma(2,2);
  }
  
  inline Box operator*(Box const &ba, double d)
  {
    Box b(ba);
    b *= d;
    return b;
  }
  
  inline Box operator+(Box const &ba, Box const &bb)
  {
    Box b(ba);
    b += bb;
    return b;
  }
  
  inline Box operator-(Box const &ba, Box const &bb)
  {
    Box b(ba);
    b -= bb;
    return b;
  }
  
  /**
   * square a box
   */
  inline Box square(Box const &ba)
  {
    Box b;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        b(i)(j) = ba(i)(j) * ba(i)(j);
    return b;
  } 
  /**
   * Box to the power of x
   */
  inline Box pow(Box const &ba, const double & x)
  {
    Box b;
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        b(i)(j) = std::pow(ba(i)(j), x);
    return b;
  }
  /**
   * Heaviside step function 
   */
  inline double heaviside(double xi)
  {
    return xi>=0.0 ? 1.0 : 0.0;
  }
} // math

#endif

