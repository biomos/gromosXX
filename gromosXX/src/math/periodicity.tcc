/**
 * @file periodicity.tcc
 * implementation of the periodic boundary condition functions.
 */

template<math::boundary_enum b>
math::Periodicity<b>::Periodicity(Matrix &box,
				  boundary_enum boundary)
  : m_box(box),
    m_boundary(boundary)
{
}

math::Periodicity<math::vacuum>::Periodicity(Matrix &box, 
					     boundary_enum boundary)
  : m_box(box)
{
  assert(boundary == vacuum); 
}

math::Periodicity<math::triclinic>::Periodicity(Matrix &box,
						boundary_enum boundary)
  : m_box(box)
{
  assert(boundary == triclinic);
}

template<math::boundary_enum b>
void math::Periodicity<b>::nearest_image(Vec const &v1, Vec const &v2,
					 Vec &nim)const
{
  switch(m_boundary){
    case vacuum:
      nim = v1 - v2;
      break;
    case triclinic:
      for(int d=0; d<3; ++d){
	nim(d) = v1(d) - v2(d);
	if (fabs(nim(d)) >= m_box(d)(d) * 0.5)
	  nim(d) -= m_box(d)(d) * rint(nim(d) / m_box(d)(d));
      }
      break;
    default:
      /*      io::messages.add("undefined boundary condition",
		       "math::Periodicity",
		       io::message::critical);
      */
      throw std::runtime_error("undefined boundary condition");
  }
}      

void math::Periodicity<math::vacuum>::nearest_image(Vec const &v1, Vec const &v2,
					      Vec &nim)const
{
  nim = v1 - v2;
}

void math::Periodicity<math::triclinic>::nearest_image(Vec const &v1, Vec const &v2,
						 Vec &nim)const
{
  for(int d=0; d<3; ++d){
    nim(d) = v1(d) - v2(d);
    if (fabs(nim(d)) >= m_box(d)(d) * 0.5)
      nim(d) -= m_box(d)(d) * rint(nim(d) / m_box(d)(d));
  }
}

template<math::boundary_enum b>
void math::Periodicity<b>::boundary(boundary_enum const boundary)
{
  m_boundary = boundary;
}

void math::Periodicity<math::vacuum>::boundary(boundary_enum const b)
{
  assert(b == vacuum);
}

void math::Periodicity<math::triclinic>::boundary(boundary_enum const b)
{
  assert(b == triclinic);
}

template<math::boundary_enum b>
math::boundary_enum const math::Periodicity<b>::boundary()const
{
  return m_boundary;
}

math::boundary_enum const math::Periodicity<math::vacuum>::boundary()const
{
  return vacuum;
}

math::boundary_enum const math::Periodicity<math::triclinic>::boundary()const
{
  return triclinic;
}
