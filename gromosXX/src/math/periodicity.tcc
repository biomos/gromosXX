/**
 * @file periodicity.tcc
 * implementation of the periodic boundary condition functions.
 */

math::periodicity<math::vacuum>::periodicity(Matrix &box)
  : m_box(box)
{
}

math::periodicity<math::triclinic>::periodicity(Matrix &box)
  : m_box(box)
{
}


void math::periodicity<math::vacuum>::nearest_image(Vec const &v1, Vec const &v2,
					      Vec &nim)
{
  nim = v1 - v2;
}

void math::periodicity<math::triclinic>::nearest_image(Vec const &v1, Vec const &v2,
						 Vec &nim)
{
  for(int d=0; d<3; ++d){
    nim(d) = v1(d) - v2(d);
    if (fabs(nim(d)) >= m_box(d)(d) * 0.5)
      nim(d) -= m_box(d)(d) * rint(nim(d) / m_box(d)(d));
  }
}
