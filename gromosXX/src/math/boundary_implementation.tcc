/**
 * @file boundary_implementation.tcc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math

#include "../debug.h"

template<math::boundary_enum b>
inline math::Boundary_Implementation<b>
::Boundary_Implementation(boundary_enum boundary)
  : m_boundary(boundary)
{
}

inline math::Boundary_Implementation<math::vacuum>
::Boundary_Implementation(boundary_enum boundary)
{
  assert(boundary == vacuum); 
}

inline math::Boundary_Implementation<math::triclinic>
::Boundary_Implementation(boundary_enum boundary)
{
  assert(boundary == triclinic);
}

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
::nearest_image(Vec const &v1, Vec const &v2,
		Vec &nim)const
{
  switch(m_boundary){
    case vacuum:
      nim = v1 - v2;
      break;
    case triclinic:
      nim = v1 - v2;

      for(int d=0; d<3; ++d){
	// i think the if statement might be wrong for really 
	// triclinic cases!
	if (fabs(nim(d)) >= m_box(d)(d) * 0.5){
	  // have to change all three components!
	  nim += m_box(d) * rint(dot(m_cross_K_L_M(d), nim));
	}
      }
      break;
    default:
      io::messages.add("undefined boundary condition",
		       "math::Boundary_Implementation",
		       io::message::critical);

      throw std::runtime_error("undefined boundary condition");
  }
}      

inline void math::Boundary_Implementation<math::vacuum>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  nim = v1 - v2;
}

inline void math::Boundary_Implementation<math::triclinic>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  nim = v1 - v2;
  for(int d=0; d<3; ++d){
    if (fabs(nim(d)) >= m_box(d)(d) * 0.5)
      nim += m_box(d) * rint(dot(m_cross_K_L_M(d), nim));
  }
}

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
::box_components(Vec const &v, Vec & n)const
{
  switch(m_boundary){
    case vacuum:
      n = 0;
      break;
    case triclinic:
      for(int d=0; d<3; ++d){
	// have to change all three components!
	n(d) = -dot(m_cross_K_L_M(d), v);
      }
      break;
    default:
      io::messages.add("undefined boundary condition",
		       "math::Boundary_Implementation",
		       io::message::critical);

      throw std::runtime_error("undefined boundary condition");
  }
}      

inline void math::Boundary_Implementation<math::vacuum>
::box_components(Vec const &v, Vec & n)const
{
  n = 0;
}

inline void math::Boundary_Implementation<math::triclinic>
::box_components(Vec const &v, Vec & n)const
{
  for(int d=0; d<3; ++d){
    n = -dot(m_cross_K_L_M(d), v);
  }
}

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
::recalc_shift_vectors(size_t const num_cells[3])
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

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
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

inline void math::Boundary_Implementation<math::triclinic>
::recalc_shift_vectors(size_t const num_cells[3])
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

inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors(size_t const num_cells[3])
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

inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].pos = 
	  0 * box()(0) +
	  0 * box()(1) +
	  0 * box()(2);
	
      }
    }
  }  
}

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
::boundary_condition(boundary_enum const boundary)
{
  m_boundary = boundary;
  if (m_boundary == triclinic){
    assert(m_volume != 0);
    m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -m_volume;
    m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / m_volume;
    m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -m_volume;
  }
  
}

inline void math::Boundary_Implementation<math::vacuum>
::boundary_condition(boundary_enum const b)
{
  assert(b == vacuum);
}

inline void math::Boundary_Implementation<math::triclinic>::
boundary_condition(boundary_enum const b)
{
  assert(b == triclinic);
}

template<math::boundary_enum b>
inline math::boundary_enum const math::Boundary_Implementation<b>
::boundary_condition()const
{
  return m_boundary;
}

inline math::boundary_enum const math::Boundary_Implementation<math::vacuum>
::boundary_condition()const
{
  return vacuum;
}

inline 
math::boundary_enum const math::Boundary_Implementation<math::triclinic>
::boundary_condition()const
{
  return triclinic;
}

// the box stuff
// -------------

// accessors
template<math::boundary_enum b>
inline math::Box const & math::Boundary_Implementation<b>::box()const
{
  return m_box;
}

inline math::Box const math::Boundary_Implementation<math::vacuum>::box()const
{
  return m_box;
}

inline math::Box const & math::Boundary_Implementation<math::triclinic>::box()const
{
  return m_box;
}

template<math::boundary_enum b>
inline double math::Boundary_Implementation<b>::volume()const
{
  return m_volume;
}

inline double math::Boundary_Implementation<math::vacuum>::volume()const
{
  return 0;
}

inline double math::Boundary_Implementation<math::triclinic>::volume()const
{
  return m_volume;
}

template<math::boundary_enum b>
inline const double math::Boundary_Implementation<b>
::box(size_t const d1, size_t const d2)const
{
  return m_box(d1)(d2);
}

inline const double math::Boundary_Implementation<math::vacuum>
::box(size_t const d1, size_t const d2)const
{
  return m_box(d1)(d2);
}

inline const double math::Boundary_Implementation<math::triclinic>
::box(size_t const d1, size_t const d2)const
{
  return m_box(d1)(d2);
}

// set the box

template<math::boundary_enum b>
inline void math::Boundary_Implementation<b>
::box(math::Box const &m)
{
  m_box = m;
  m_volume = dot(cross(m_box(K), m_box(L)), m_box(M));
  // std::cout << "m_volume = " << m_volume << std::endl;
  
  if (m_boundary==triclinic){
    assert(m_volume != 0);
    m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -m_volume;
    m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / m_volume;
    m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -m_volume;
  }
}

inline void math::Boundary_Implementation<math::vacuum>
::box(math::Box const &m)
{
  m_box = m;
}

inline void math::Boundary_Implementation<math::triclinic>
::box(math::Box const &m)
{
  m_box = m;
  m_volume = dot(cross(m_box(K), m_box(L)), m_box(M));
  
  m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -m_volume;
  m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / m_volume;
  m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -m_volume;

}

template<math::boundary_enum b>
inline typename math::Boundary_Implementation<b>::shift_struct &
math::Boundary_Implementation<b>
::shift(size_t const i)
{
  assert(27 > i);
  return m_shift[i];
}

template<math::boundary_enum b>
inline typename math::Boundary_Implementation<b>::shift_struct const &
math::Boundary_Implementation<b>
::shift(size_t const i)const
{
  assert(27 > i);
  return m_shift[i];
}

inline math::Boundary_Implementation<math::vacuum>::shift_struct &
math::Boundary_Implementation<math::vacuum>
::shift(size_t const i)
{
  assert(27 > i);
  return m_shift[i];
}

inline math::Boundary_Implementation<math::vacuum>::shift_struct const &
math::Boundary_Implementation<math::vacuum>
::shift(size_t const i)const
{
  assert(27 > i);
  return m_shift[i];
}

inline math::Boundary_Implementation<math::triclinic>::shift_struct &
math::Boundary_Implementation<math::triclinic>
::shift(size_t const i)
{
  assert(27 > i);
  return m_shift[i];
}

inline math::Boundary_Implementation<math::triclinic>::shift_struct const &
math::Boundary_Implementation<math::triclinic>
::shift(size_t const i)const
{
  assert(27 > i);
  return m_shift[i];
}

