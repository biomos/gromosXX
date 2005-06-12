/**
 * @file gmath.h
 * mathematical definitions.
 */

#ifndef INCLUDED_MATH_H
#define INCLUDED_MATH_H

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
  };
  
#endif

