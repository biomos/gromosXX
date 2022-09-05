/**
 * @file  transformation.cc
 * implementation of box transformations
 */

#include "../stdheader.h"
#include "../math/transformation.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

/**
 * calculate the rotation matrix R from the box
 */
math::Matrixl math::rmat(math::Box const & box) {
  const math::Vecl Rx = box(0) / math::abs(box(0));
  /* is this wrong???
   * const math::Vecl Ry_aux = box(1)
          - math::dot(box(0), box(1)) * box(0)
          / (math::abs(box(0)) * math::abs(box(1)));*/
  const math::Vecl Ry_aux = box(1)
          - math::dot(box(0), box(1)) * box(0)
          / (math::abs(box(0))*math::abs(box(0)));
  const math::Vecl Ry = Ry_aux / math::abs(Ry_aux);
  const math::Vecl Rz = math::cross(Rx, Ry);
  math::Matrixl R(Rx, Ry, Rz);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (std::abs(R(i, j)) <= math::epsilon)
        R(i, j) = 0.0;
  return R;
}
math::Matrix math::rmat(math::Matrix const & box) {
  math::Vec a, b, c;
  for(int i=0; i<3; i++){
    a(i)=box(0,i);
    b(i)=box(1,i);
    c(i)=box(2,i);
  }
        
   const math::Vecl Rx = a / math::abs(a);
  /* is this wrong???
   * const math::Vecl Ry_aux = box(1)
          - math::dot(box(0), box(1)) * box(0)
          / (math::abs(box(0)) * math::abs(box(1)));*/
  const math::Vecl Ry_aux = b
          - math::dot(a, b) * a
          / math::abs2(a);
  const math::Vecl Ry = Ry_aux / math::abs(Ry_aux);
  const math::Vecl Rz = math::cross(Rx, Ry);
  math::Matrix R(Rx, Ry, Rz, true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (fabs(R(i, j)) <= math::epsilon)
        R(i, j) = 0.0;
  return R;
}
/**
 * calculate the rotation matrix R from the angles
 */
math::Matrixl math::rmat(double const & phi, double const & theta,
        double const & psi) {
  const math::Vecl Rx(cosl(theta) * cosl(phi),
          cosl(theta) * sinl(phi),
          -sinl(theta));
  const math::Vecl Ry(sinl(psi) * sinl(theta) * cosl(phi) - cosl(psi) * sinl(phi),
          sinl(psi) * sinl(theta) * sinl(phi) + cosl(psi) * cosl(phi),
          sinl(psi) * cosl(theta));
  const math::Vecl Rz(cosl(psi) * sinl(theta) * cosl(phi) + sinl(psi) * sinl(phi),
          cosl(psi) * sinl(theta) * sinl(phi)+(-sinl(psi) * cosl(phi)),
          cosl(psi) * cosl(theta));

  //math::Matrixl R(Rx, Ry, Rz, true);
  math::Matrixl R(Rx, Ry, Rz);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (std::abs(R(i, j)) <= math::epsilon)
        R(i, j) = 0.0;
  return R;
}

/**
 * calculate the transformation matrix S.
 */
math::Matrixl math::smat(math::Box const & box, math::boundary_enum const b) {
  if (b == math::rectangular) {
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(0.0, 1.0, 0.0);
    math::Vecl Sz(0.0, 0.0, 1.0);
    math::Matrixl S(Sx, Sy, Sz);
    return S;
  } else {
    long double alpha = 0.0, beta = 0.0, gamma = 0.0;
    alpha = acos(math::costest(dot(box(1), box(2)) / (abs(box(1)) * abs(box(2)))));
    beta = acos(math::costest(dot(box(0), box(2)) / (abs(box(0)) * abs(box(2)))));
    gamma = acos(math::costest(dot(box(0), box(1)) / (abs(box(0)) * abs(box(1)))));

    long double cosdelta = (cosl(alpha) - cosl(beta) * cosl(gamma)) / (sinl(beta) * sinl(gamma));
    long double sindelta = sqrtl(1 - cosdelta * cosdelta);
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(cosl(gamma),
            sinl(gamma),
            0.0);
    math::Vecl Sz(cosl(beta),
            cosdelta * sinl(beta),
            sindelta * sinl(beta));
    math::Matrixl S(Sx, Sy, Sz);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (std::abs(S(i, j)) < math::epsilon)
          S(i, j) = 0.0;
    return S;
  }
}

/**
 * calculate the inverse of the transformation matrix S.
 */
math::Matrixl math::sinvmat(math::Box const & box, math::boundary_enum const b) {
  if (b == math::rectangular) {
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(0.0, 1.0, 0.0);
    math::Vecl Sz(0.0, 0.0, 1.0);
    math::Matrixl Sinv(Sx, Sy, Sz);
    return Sinv;
  } else {
    long double alpha = 0.0, beta = 0.0, gamma = 0.0;
    alpha = acos(math::costest(dot(box(1), box(2)) / (abs(box(1)) * abs(box(2)))));
    beta = acos(math::costest(dot(box(0), box(2)) / (abs(box(0)) * abs(box(2)))));
    gamma = acos(math::costest(dot(box(0), box(1)) / (abs(box(0)) * abs(box(1)))));

    long double cosdelta = (cosl(alpha) - cosl(beta) * cosl(gamma)) / (sinl(beta) * sinl(gamma));
    long double sindelta = sqrtl(1 - cosdelta * cosdelta);
    long double cotdelta = cosdelta / sindelta;
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(-cotdelta,
            1 / sinl(gamma), 0.0);
    math::Vecl Sz((cotdelta / tanl(gamma) - 1 / (sindelta * tanl(beta))),
            -cotdelta / sinl(gamma),
            1 / (sinl(beta) * sindelta));
    math::Matrixl Sinv(Sx, Sy, Sz);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (Sinv(i, j) == -0.0)
          Sinv(i, j) = 0.0;
    return Sinv;
  }
}

/**
 * calculate the transformation matrix M.
 */
math::Matrixl math::mmat(math::Box const & box, math::boundary_enum const b) {
  return math::product(math::rmat(box), math::smat(box, b));
}

/**
 * calculate the inverse of the transformation matrix M.
 */
math::Matrixl math::minvmat(math::Box const & box, math::boundary_enum const b) {
  return math::product(math::sinvmat(box, b), math::transpose(math::rmat(box)));
}

math::Matrix math::truncoct_triclinic_rotmat(bool forward) {
  static const double sq3 = sqrt(3.0);
  static const double sq3i = 1.0/sq3;
  static const double sq2i = 1.0/sqrt(2.0);

  math::Matrix rot(Vec(sq3i, -2*sq2i*sq3i, 0),
	       Vec(sq3i, sq3i*sq2i, -sq2i),
	       Vec(sq3i, sq2i*sq3i, sq2i));

  if (!forward)
    rot = math::transpose(rot);

  return rot;
}

void math::truncoct_triclinic_box(math::Box & box, bool forward) {
  static const double third = 1.0 / 3.0;
  static const double sq3 = sqrt(3.0);
  static const double sq3i = 1.0/sq3;

  if(forward){
    const double d = 0.5 *sq3 * box(0)(0);
    const double sq2 = sqrt(2.0);

    box(0) = math::Vec(d, 0.0, 0.0);
    box(1) = math::Vec(third * d, 2.0 * third * sq2 * d, 0.0);
    box(2) = math::Vec(-third * d, third * sq2 * d, third * sqrt(6.0) * d);
  }  else {
    const double d = 2.0 * box(0)(0) * sq3i;
    box(0) = math::Vec(d, 0.0, 0.0);
    box(1) = math::Vec(0.0, d, 0.0);
    box(2) = math::Vec(0.0, 0.0, d);
  }
}

void math::truncoct_triclinic(math::VArray & pos, bool forward) {
  DEBUG(6, (forward ? "truncated to triclinic" : "triclinic to truncated"));
  math::Matrix rot = math::truncoct_triclinic_rotmat(forward);
  DEBUG(8, "rotation matrix:\n\t" << math::m2s(rot));

  // apply transformation
  for(unsigned int i = 0; i < pos.size(); ++i) {
    pos(i) = math::product(rot, pos(i));
  }
}
