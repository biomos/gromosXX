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
      // std::cout << "nim: " << nim << std::endl;
      // std::cout << "m_box: " << m_box << std::endl;

      for(int d=0; d<3; ++d){
	if (fabs(nim(d)) >= m_box(d)(d) * 0.5){
	  // have to change all three components!
	  nim += m_box(d) * rint(dot(m_cross_K_L_M(d), nim));
	}
	// std::cout << "nim (" << d << "): " << nim << std::endl;	
      }
      break;
    default:
      /*
      io::messages.add("undefined boundary condition",
		       "math::Boundary_Implementation",
		       io::message::critical);
      */
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
inline math::Matrix const & math::Boundary_Implementation<b>::box()const
{
  return m_box;
}

inline math::Matrix const math::Boundary_Implementation<math::vacuum>::box()const
{
  return m_box;
}

inline math::Matrix const & math::Boundary_Implementation<math::triclinic>::box()const
{
  return m_box;
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
::box(math::Matrix const &m)
{
  m_box = m;
  m_volume = dot(cross(m_box(K), m_box(L)), m_box(M));

  if (m_boundary==triclinic){
    assert(m_volume != 0);
    m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -m_volume;
    m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / m_volume;
    m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -m_volume;
  }
}

inline void math::Boundary_Implementation<math::vacuum>
::box(math::Matrix const &m)
{
  m_box = m;
}

inline void math::Boundary_Implementation<math::triclinic>
::box(math::Matrix const &m)
{
  m_box = m;
  m_volume = dot(cross(m_box(K), m_box(L)), m_box(M));
  
  m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -m_volume;
  m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / m_volume;
  m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -m_volume;

}

  
  
