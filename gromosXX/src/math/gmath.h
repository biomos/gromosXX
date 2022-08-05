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
    GenericVec() : GenericVec(0) {}
    explicit GenericVec(numeric_type d) { d_v[0] = d_v[1] = d_v[2] = d; }
    template<typename numeric_type_b>
    GenericVec(GenericVec<numeric_type_b> const & v) { d_v[0] = v(0); d_v[1] = v(1); d_v[2] = v(2); }
    GenericVec(numeric_type d1, numeric_type d2, numeric_type d3) { d_v[0] = d1; d_v[1] = d2; d_v[2] = d3; }

    numeric_type operator()(int i)const { assert(i>=0 && i<3); return d_v[i]; }
    numeric_type & operator()(int i) { assert(i >= 0 && i < 3); return d_v[i]; }
    numeric_type & operator[](unsigned int i) { assert(i < 3); return d_v[i]; }
    const numeric_type & operator[](unsigned int i) const { assert(i < 3); return d_v[i]; }
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
    
    inline GenericVec<numeric_type> norm() const;
    
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

    /**
     * Subtract from every Vector in VArray a vector
     */
    VArray & operator-=(Vec const &v)
    {
      VArray::iterator it = this->begin(),
                       to = this->end();
      for (; it != to; it++){
        *it -= v;
      }
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
    SArray & operator=(double d)
    {
      for(std::vector<double>::iterator it=begin(), to=end(); it!=to; ++it)
	*it = d;
      return *this;
    }
  };
  
  // the box should be known to the Matrix
  class Box;
  
  
  template<typename numeric_type> class GenericSymmetricMatrix;
  /**
   * Matrix.
   */
  template<typename numeric_type>
  class GenericMatrix
  {
  private:
    numeric_type m[3][3];
  public:
    GenericMatrix() : GenericMatrix(0) {}
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
    inline GenericMatrix(const GenericSymmetricMatrix<numeric_type_b> & ma);

    template<typename numeric_type_b>
    inline GenericMatrix(
            numeric_type_b d1, numeric_type_b d2, numeric_type_b d3,
            numeric_type_b d4, numeric_type_b d5, numeric_type_b d6,
            numeric_type_b d7, numeric_type_b d8, numeric_type_b d9) {
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
    GenericMatrix & operator=(const GenericMatrix<numeric_type_b> & mat)
    {
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  m[i][j] = mat(i,j);
      return *this;
    }

    inline
    numeric_type operator()(int i, int j)const 
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      return m[i][j];
    }

    inline
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
    inline GenericMatrix<numeric_type> & operator+=(GenericSymmetricMatrix<numeric_type_b> const & mat);

    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator-=(GenericMatrix<numeric_type_b> const & mat) {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          m[i][j] -= mat(i, j);
      return *this;
    }
    
    template<typename numeric_type_b>
    inline GenericMatrix<numeric_type> & operator-=(GenericSymmetricMatrix<numeric_type_b> const & mat);
    
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
    
    inline Vec operator*(const Vec v){
      Vec product = Vec(0);
      for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
          product[i] += m[i][j] * v[j];
        }
      }
      return product;
    }
  };

  /**
   * @class GenericSymmetricMatrix
   * @ingroup math
   *
   * Matrix, a generic and symmetric version
   * it looks like this:
   * @verbatim
  [[  m[0]  m[1]  m[2]  ]
   [  m[1]  m[3]  m[4]  ]
   [  m[2]  m[4]  m[5]  ]]
     @endverbatim
   */
  template<typename numeric_type>
  class GenericSymmetricMatrix
  {
  private:
    /**
     * the 6 values of the symmetric matrix. See @ref math::GenericSymmetricMatrix
     * for details about the elements.
     */
    numeric_type m[6];
  public:
    /**
     * constructor
     */
    GenericSymmetricMatrix() : GenericSymmetricMatrix(0) {}
    
    /**
     * construct it from a numerical value
     */
    template<typename numeric_type_b>
    explicit GenericSymmetricMatrix(numeric_type_b d) 
    {
      m[0] = m[1] = m[2] = m[3] = m[4] = m[5] = d;
    }

    /**
     * construct it from 6 numerical values
     */
    template<typename numeric_type_b>
    inline GenericSymmetricMatrix(
            const numeric_type_b &d1, const numeric_type_b &d2, const numeric_type_b &d3,
            const numeric_type_b &d4, const numeric_type_b &d5, const numeric_type_b &d6) {
      m[0] = d1;
      m[1] = d2;
      m[2] = d3;
      m[3] = d4;
      m[4] = d5;
      m[5] = d6;
    }
   
    /**
     * assignment operator for a single numerical value
     */
    template<typename numeric_type_b>
    GenericSymmetricMatrix & operator=(numeric_type_b d)
    {
      m[0] = m[1] = m[2] = m[3] = m[4] = m[5] = d;
      return *this;
    }
    
    /**
     * assignment operator for a @ref math::GenericMatrix
     */
    template<typename numeric_type_b>
    GenericSymmetricMatrix & operator=(const GenericMatrix<numeric_type_b> & mat)
    {
      m[0] = mat(0,0);
      m[1] = mat(0,1);
      m[2] = mat(0,2);
      m[3] = mat(1,1);
      m[4] = mat(1,2);
      m[5] = mat(2,2);
      return *this;
    }
    
    /**
     * this operator is used to access a matrix element.
     * It is slow and should be avoided. const version
     */
    numeric_type operator()(int i, int j)const 
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      if (j < i) { // swap the indices
        int tmp = i;
        i = j; j = tmp;
      }
      if (i == 0)
        return m[j];
      else if (i == 1)
        return m[2+j];
      else
        return m[5];
    }
    
    /**
     * this operator is used to access a matrix element.
     * It is slow and should be avoided. 
     */
    numeric_type & operator()(int i, int j)
    {
      assert( i>=0 && i<3 && j>=0 && j<3 );
      if (j < i) { // swap the indices
        int tmp = i;
        i = j; j = tmp;
      }
      if (i == 0)
        return m[j];
      else if (i == 1)
        return m[2+j];
      else
        return m[5];
    }
    
    /**
     * access one of the six matrix elements
     * const version
     */
    numeric_type operator()(int i) const {
      assert(i >= 0 && i < 6);
      return m[i];
    }
    
    /**
     * access one of the six matrix elements
     */
    numeric_type & operator()(int i) {
      assert(i >= 0 && i < 6);
      return m[i];
    }

    /**
     * add another symetric matrix
     */
    template<typename numeric_type_b>
    inline GenericSymmetricMatrix<numeric_type> & operator+=(GenericSymmetricMatrix<numeric_type_b> const & mat) {
      for (int i = 0; i < 6; ++i)
        m[i] += mat(i);
      return *this;
    }

    /**
     * subtract another symmetric matrix
     */
    template<typename numeric_type_b>
    inline GenericSymmetricMatrix<numeric_type> & operator-=(GenericSymmetricMatrix<numeric_type_b> const & mat) {
      for (int i = 0; i < 6; ++i)
        m[i] -= mat(i);
      return *this;
    }
    
    /**
     * scale the matrix
     */
    template<typename numeric_type_b>
    inline GenericSymmetricMatrix<numeric_type> & operator*=(const numeric_type_b & d) {
      for (int i = 0; i < 6; ++i)
        m[i] *= d;
      return *this;
    }
    
    /**
     * scale the matrix (by division)
     */
    template<typename numeric_type_b>
    inline GenericSymmetricMatrix<numeric_type> & operator/=(const numeric_type_b & d) {
      for (int i = 0; i < 6; ++i)
        m[i] /= d;
      return *this;
    }
    
    /**
     * add elements to the diagonal of the matrix
     */
    template<typename numeric_type_b>
    inline void add_to_diagonal(const numeric_type_b & d) {
      m[0] += d; m[3] += d; m[5] += d;
    }
  };  
  
  /**
   * convert the symmetric matrix into a normal matrix
   */
  template<typename numeric_type>
  template<typename numeric_type_b>
  GenericMatrix<numeric_type>::GenericMatrix(const GenericSymmetricMatrix<numeric_type_b> & ma) {
    m[0][0] = ma(0);
    m[0][1] = m[1][0] = ma(1);
    m[0][2] = m[2][0] = ma(2);
    m[1][1] = ma(3);
    m[1][2] = m[2][1] = ma(4);
    m[2][2] = ma(5);   
  }
  
  /**
   * a double matrix
   */
  typedef GenericMatrix<double> Matrix;
  /**
   * a long double matrix
   */
  typedef GenericMatrix<long double> Matrixl;

  /**
   * a double symmetric matrix
   */
  typedef GenericSymmetricMatrix<double> SymmetricMatrix;
  
  /**
   * add a symmetric matrix to a normal matrix - slow
   */
  template<typename numeric_type>
  template<typename numeric_type_b>
  inline GenericMatrix<numeric_type> & GenericMatrix<numeric_type>::operator+=(GenericSymmetricMatrix<numeric_type_b> const & mat) {
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        m[i][j] += mat(i ,j);

    return *this;
  }
  
  /**
   * subtract a symmetric matrix from a normal matrix - slow
   */
  template<typename numeric_type>
  template<typename numeric_type_b>
  inline GenericMatrix<numeric_type> & GenericMatrix<numeric_type>::operator-=(GenericSymmetricMatrix<numeric_type_b> const & mat) {
    for(int i = 0; i < 3; i++)
      for(int j = i; j < 3; j++)
        m[i][j] -= mat(i ,j);

    return *this;
  }
  
  /**
   * Box.
   */
  class Box
  {
  private:
    Vec d_b[3];
  public:
    Box() {
     // avoid occurrences of uninitialized vacuum boxes
     d_b[0] = 0.0;
     d_b[1] = 0.0;
     d_b[2] = 0.0;
    }
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
   Box(const Matrixl &m)
    {
       for (unsigned int i = 0; i < 3; ++i) {
        d_b[i]=(m(0,i),m(1,i),m(2,i));

       }   
    }
    Box(const Box& other)
    {
      // copy constructor
      d_b[0] = other.d_b[0];
      d_b[1] = other.d_b[1];
      d_b[2] = other.d_b[2];
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
    truncoct = -1, // 01.12.2009, Andreas: set from 3 to -1
    /**
     * non-specialized version
     */
    any = 4 // 01.12.2009, Andreas: could also be set to 3 but I leave it at 4
  };
  /**
   * @enum virial_enum
   * virial enum
   */
  enum virial_enum { 
    /** no virial */ 
    no_virial = 0, 
    /** molecular/group virial */
    molecular_virial = 2, 
    /** atomic virial */
    atomic_virial = 1 
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
    pcouple_full_anisotropic = 3,
    /** semianisotropic pressure scaling */
    pcouple_semi_anisotropic = 4
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
   * speed of light.
   */
  extern double spd_l;

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

  /**
   * Copies a three-dimensional vector of doubles into a one-dimensional vector with the specified offset and scaling. 
   * Useful to prepare data for Fortran API calls. No bounce checks.
   * 
   * @param v_f One-dimensional vector of doubles ready to go into the Fortran code
   * @param v_c Cartesian coordinates of an atom / gradient to be copied
   * @param offset The specified offset
   * @param scaling The scaling factor
   */
  inline void vector_c2f(std::vector<double>& v_f
                , const math::Vec& v_c
                , const unsigned int offset
                , const double scaling) 
  {
    for (unsigned int i = 0; i < 3; ++i) {
      v_f[3 * offset + i] = v_c(i) * scaling;
    }
  }

  /**
   * Copies a three-dimensional vector of doubles into a one-dimensional vector with the specified offset. 
   * Useful to prepare data for Fortran API calls. No bounce checks.
   * 
   * @param v_f One-dimensional vector of doubles ready to go into the Fortran code
   * @param v_c Cartesian coordinates of an atom / gradient to be copied
   * @param offset The specified offset
   */
  inline void vector_c2f(std::vector<double>& v_f
                , const math::Vec& v_c
                , const unsigned int offset) 
  {
    vector_c2f(v_f, v_c, offset, 1.0);
  }

  /**
   * Copies a one-dimensional vector of doubles into a three-dimensional vector with the specified offset and scaling. 
   * Useful to process data received from Fortran API calls. No bounce checks.
   * 
   * @param v_f One dimensional vector of doubles received from a Fortran API call
   * @param v_c Cartesian coordinates of an atom / gradient to be copied to
   * @param offset The specified offset
   * @param scaling The scaling factor
   */
  inline void vector_f2c(const std::vector<double>& v_f
                , math::Vec& v_c
                , const unsigned int offset
                , const double scaling) 
  {
    for (unsigned int i = 0; i < 3; ++i) {
      v_c(i) = v_f[3 * offset + i] * scaling;
    }
  }

  /**
   * Copies a one-dimensional vector of doubles into a three-dimensional vector with the specified offset. 
   * Useful to process data received from Fortran API calls. No bounce checks.
   * 
   * @param v_f One dimensional vector of doubles received from a Fortran API call
   * @param v_c Cartesian coordinates of an atom / gradient to be copied to
   * @param offset The specified offset
   */
  inline void vector_f2c(const std::vector<double>& v_f
                , math::Vec& v_c
                , const unsigned int offset) 
  {
    vector_f2c(v_f, v_c, offset, 1.0);
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
  
  /**
   * calculate the outer product of two vectors
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_b> tensor_product(GenericVec<numeric_type_b> const &v1, GenericVec<numeric_type_a> const &v2)
  {
    GenericMatrix<numeric_type_b> result;
    for (unsigned int a=0; a<3; ++a){            
       for(unsigned int b=0; b<3; ++b)
         result(b, a) = v1(b) * v2(a);
    }
    return result;
  }
  
  /**
   * calculate the ouer product of two vectors
   *
   * This version asumes that the resulting matrix is symmetrical without
   * checking whether this is true!
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_b> symmetric_tensor_product(GenericVec<numeric_type_b> const &v1, GenericVec<numeric_type_a> const &v2)
  {
    GenericSymmetricMatrix<numeric_type_b> result;
    result(0) = v1(0) * v2(0);
    result(1) = v1(0) * v2(1);
    result(2) = v1(0) * v2(2);
    result(3) = v1(1) * v2(1);
    result(4) = v1(1) * v2(2);
    result(5) = v1(2) * v2(2);
    return result;
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
  template<typename numeric_type>
  inline Box product(GenericMatrix<numeric_type> const &m1, Box const &m2)
  {
    Box m(0.0);
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	for(int k=0; k<3; ++k)
	  m(i)(j) += m1(i,k) * m2(k)(j);
    return m;
  }
  template<typename numeric_type>
  inline Matrix product(Box const &m1, GenericMatrix<numeric_type> const &m2)
  {
    Matrix m(0.0);   
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	for(int k=0; k<3; ++k)
	  m(i,j) += m1(i)(k) * m2(k,j);
    return m;
  }
  template<typename numeric_type>
  inline std::string m2s(GenericMatrix<numeric_type> const & m)
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
  
  inline std::string m2s(GenericSymmetricMatrix<double> const & m)
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
  inline GenericMatrix<numeric_type_a> operator*(numeric_type_b d, GenericMatrix<numeric_type_a> const &ma)
  {
    return GenericMatrix<numeric_type_a>(ma(0,0) * d, ma(0,1) * d, ma(0,2) * d,
            ma(1,0) * d, ma(1,1) * d, ma(1,2) * d,
            ma(2,0) * d, ma(2,1) * d, ma(2,2) * d);
  }  
  
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericMatrix<numeric_type_a> operator/(GenericMatrix<numeric_type_a> const &ma, numeric_type_b d)
  {
    return GenericMatrix<numeric_type_a>(ma(0,0) / d, ma(0,1) / d, ma(0,2) / d,
            ma(1,0) / d, ma(1,1) / d, ma(1,2) / d,
            ma(2,0) / d, ma(2,1) / d, ma(2,2) / d);
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
   * scale a symmetric matrix
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_a> operator*(GenericSymmetricMatrix<numeric_type_a> const &ma, numeric_type_b d)
  {
    return GenericSymmetricMatrix<numeric_type_a>(ma(0) * d, ma(1) * d, ma(2) * d,
            ma(3) * d, ma(4) * d, ma(5) * d);
  }
  
  /**
   * scale a symmetric matrix
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_a> operator*(numeric_type_b d, GenericSymmetricMatrix<numeric_type_a> const &ma)
  {
    return GenericSymmetricMatrix<numeric_type_a>(ma(0) * d, ma(1) * d, ma(2) * d,
            ma(3) * d, ma(4) * d, ma(5) * d);
  }  
  
  /**
   * scale a symmetric matrix
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_a> operator/(GenericSymmetricMatrix<numeric_type_a> const &ma, numeric_type_b d)
  {
    return GenericSymmetricMatrix<numeric_type_a>(ma(0) / d, ma(1) / d, ma(2) / d,
            ma(3) / d, ma(4) / d, ma(5) / d);
  }  
  
  /**
   * add two symmetric matrices
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_a> operator+(GenericSymmetricMatrix<numeric_type_a> const &ma, GenericSymmetricMatrix<numeric_type_b> const &mb)
  {
    GenericSymmetricMatrix<numeric_type_a> m;
    for(int i = 0; i < 6; ++i)
      m(i) = ma(i) + mb(i);
    return m;
  }
  
  /**
   * subtract two symmetric matrices
   */
  template<typename numeric_type_a, typename numeric_type_b>
  inline GenericSymmetricMatrix<numeric_type_a> operator-(GenericSymmetricMatrix<numeric_type_a> const &ma, GenericSymmetricMatrix<numeric_type_b> const &mb)
  {
    GenericSymmetricMatrix<numeric_type_a> m;
    for(int i = 0; i < 6; ++i)
      m(i) = ma(i) - mb(i);
    return m;
  }
  
  /**
   * square a matrix
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
    const double det_ma = 1.0 / det(ma);
    return GenericMatrix<numeric_type > (
            (ma(1, 1) * ma(2, 2) - ma(1, 2) * ma(2, 1)) * det_ma,
            (ma(0, 2) * ma(2, 1) - ma(0, 1) * ma(2, 2)) * det_ma,
            (ma(0, 1) * ma(1, 2) - ma(0, 2) * ma(1, 1)) * det_ma,
            (ma(1, 2) * ma(2, 0) - ma(1, 0) * ma(2, 2)) * det_ma,
            (ma(0, 0) * ma(2, 2) - ma(0, 2) * ma(2, 0)) * det_ma,
            (ma(0, 2) * ma(1, 0) - ma(0, 0) * ma(1, 2)) * det_ma,
            (ma(1, 0) * ma(2, 1) - ma(1, 1) * ma(2, 0)) * det_ma,
            (ma(0, 1) * ma(2, 0) - ma(0, 0) * ma(2, 1)) * det_ma,
            (ma(0, 0) * ma(1, 1) - ma(0, 1) * ma(1, 0)) * det_ma);
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
    return ma(0,0) + ma(1,1) + ma(2,2);
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
  
  /**
   * norm of a vector
   */
  template<typename numeric_type>
  GenericVec<numeric_type> GenericVec<numeric_type>::norm() const {
    return *this / abs(*this);
  }
  
  /**
   * isnan function 
   */
  inline bool gisnan(double x) {
      return (x) != (x);
  }
  
} // math

#endif /* INCLUDED_MATH_H */