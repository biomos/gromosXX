/**
 * @file periodicity.tcc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math

#include "../debug.h"

template<math::boundary_enum b>
math::Periodicity<b>::Periodicity(boundary_enum boundary)
  : Boundary_Implementation<b>(boundary)
{
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_box(math::Vec &v)const
{
  Vec o(0, 0, 0);
  nearest_image(v, o, v);
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_positive_box(math::Vec &v)const
{
  Vec o(m_box(0)(0), m_box(1)(1), m_box(2)(2));
  o /= 2;
  nearest_image(v, o, v);
  v += o;  
}
