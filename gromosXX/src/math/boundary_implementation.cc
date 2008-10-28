/**
 * @file boundary_implementation.cc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math

#include "../util/debug.h"

#ifdef WIN32
// Converts a floating point value to an integer, very fast.
inline int rint(float param)
{
	// Uses the FloatToInt functionality
	int a;
	int *int_pointer = &a;

	__asm  fld  param
	__asm  mov  edx,int_pointer
	__asm  FRNDINT
	__asm  fistp dword ptr [edx];

	return a;
}
#endif

/*
 * Constructor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>
::Boundary_Implementation(math::Box const & b)
  : m_box(b)
{
}

/**
 * Constructor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>
::Boundary_Implementation(math::Box const & b)
 : m_box(b)
{
  for(int i=0; i<3; ++i)
    m_half_box(i) = 0.5 * abs(m_box(i));
}

/**
 * Constructor : triclinic
 */
inline math::Boundary_Implementation<math::triclinic>
::Boundary_Implementation(math::Box const & b)
  : m_box(b)
{
  double volume = dot(cross(m_box(K), m_box(L)), m_box(M));
  
  assert(volume != 0);
  m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -volume;
  m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / volume;
  m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -volume;
}

/**
 * Constructor : truncated octahedron
 */
inline math::Boundary_Implementation<math::truncoct>
::Boundary_Implementation(math::Box const & b)
 : m_box(b)
{
  for(int i=0; i<3; ++i)
    m_half_box(i) = 0.5 * abs(m_box(i));
}

////////////////////////////////////////////////////////////////////////////////
// box / shift vector accessors
////////////////////////////////////////////////////////////////////////////////

/**
 * const box accessor : vacuum
 */
inline math::Box const math::Boundary_Implementation<math::vacuum>::box()const
{
  return m_box;
}

/**
 * const box accessor : rectangular
 */
inline math::Box const & math::Boundary_Implementation<math::rectangular>::box()const
{
  return m_box;
}

/**
 * const box accessor : triclinic
 */
inline math::Box const & math::Boundary_Implementation<math::triclinic>::box()const
{
  return m_box;
}

/**
 * const box accessor : truncated octahedron
 */
inline math::Box const & math::Boundary_Implementation<math::truncoct>::box()const
{
  return m_box;
}

/**
 * box element accessor (d1,d2) : vacuum
 */
inline double math::Boundary_Implementation<math::vacuum>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * box element accessor (d1,d2) : rectangular
 */
inline double math::Boundary_Implementation<math::rectangular>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * box element accessor (d1,d2) : triclinic
 */
inline double math::Boundary_Implementation<math::triclinic>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * box element accessor (d1,d2) : truncated octahedron
 */
inline double math::Boundary_Implementation<math::truncoct>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * shift struct accessor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>::shift_struct &
math::Boundary_Implementation<math::vacuum>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>::shift_struct const &
math::Boundary_Implementation<math::vacuum>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * shift struct accessor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>::shift_struct &
math::Boundary_Implementation<math::rectangular>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>::shift_struct const &
math::Boundary_Implementation<math::rectangular>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * shift struct accessor : triclininc
 */
inline math::Boundary_Implementation<math::triclinic>::shift_struct &
math::Boundary_Implementation<math::triclinic>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : triclinic
 */
inline math::Boundary_Implementation<math::triclinic>::shift_struct const &
math::Boundary_Implementation<math::triclinic>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

////////////////////////////////////////////////////////////////////////////////
// nearest image functions
////////////////////////////////////////////////////////////////////////////////

/**
 * nearest image : vacuum
 */
inline void math::Boundary_Implementation<math::vacuum>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  nim = v1 - v2;
}

/**
 * nearest image : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::nearest_image(Vec const &v1, Vec const &v2,
		Vec &nim)const
{
  // nim = v1 - v2;

  for(int d=0; d<3; ++d){
    nim(d) = v1(d) - v2(d);

    if (fabs(nim(d)) >= m_half_box(d)){
      nim(d) -= abs(m_box(d)) * rint(nim(d)/abs(m_box(d)));

    }
  }
}      

/**
 * nearest image : triclinic
 */
inline void math::Boundary_Implementation<math::triclinic>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  // nim has to be out of the loop here!
  nim = v1 - v2;
  for(int d=0; d<3; ++d){
    // i think the if statement might be wrong for really 
    // triclinic cases!
    if (fabs(nim(d)) >= abs(m_box(d)) * 0.5)
      nim += m_box(d) * rint(dot(m_cross_K_L_M(d), nim));
  }
}

/**
 * nearest image : truncated octahedron
 */
inline void math::Boundary_Implementation<math::truncoct>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  for(int d=0; d<3; ++d){
    nim(d) = v1(d) - v2(d);

    if (fabs(nim(d)) >= m_half_box(0)){
      nim(d) -= abs(m_box(0)) * rint(nim(d)/abs(m_box(0)));

    }
  }
  if (fabs(nim(0)) + fabs(nim(1)) + fabs(nim(2)) > 0.75 *abs(m_box(0))){
    for(int d=0; d<3; ++d){
      if (nim(d) < 0)
	nim(d) += m_half_box(0);
      else
	nim(d) -= m_half_box(0);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// grid stuff
////////////////////////////////////////////////////////////////////////////////


/**
 * calculate box components of vector v : vacuum
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::vacuum>
::box_components(Vec const &v, Vec & n)const
{
  n = 0;
}

/**
 * calculate box components of vector v : rectangular
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::rectangular>
::box_components(Vec const &v, Vec & n)const
{
  for(int d=0; d<3; ++d){
    n(d) = v(d) / m_box(d)(d);
  }
}

/**
 * calculate box components of vector v : triclinic
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::triclinic>
::box_components(Vec const &v, Vec & n)const
{
  for(int d=0; d<3; ++d){
    n(d) = -dot(m_cross_K_L_M(d), v);
  }
}

/**
 * recalc shift vectors
 * and also update cell index shifts : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = k * num_cells[0];
	m_shift[index].cell[1] = l * num_cells[1];
	m_shift[index].cell[2] = m * num_cells[2];

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){
	m_shift[index].pos(0) = k * box(0,0);
	m_shift[index].pos(1) = l * box(1,1);
	m_shift[index].pos(2) = m * box(2,2);
      }
    }
  }  
}

/**
 * recalc shift vectors
 * and also update cell index shifts : triclinic
 */
inline void math::Boundary_Implementation<math::triclinic>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = k * num_cells[0];
	m_shift[index].cell[1] = l * num_cells[1];
	m_shift[index].cell[2] = m * num_cells[2];

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors : triclinic
 */
inline void math::Boundary_Implementation<math::triclinic>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors: vacuum
 * (everything is 0)
 */
inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = 0;
	m_shift[index].cell[1] = 0;
	m_shift[index].cell[2] = 0;

	m_shift[index].pos = 
	  0 * box()(0) +
	  0 * box()(1) +
	  0 * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vector : vacuum
 * (everything 0)
 */
inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].pos = 0;
	
      }
    }
  }  
}

