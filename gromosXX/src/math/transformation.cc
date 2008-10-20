/**
 * @file  transformation.cc
 * implementation of box transformations
 */

#include <stdheader.h>
#include "transformation.h"
#include "gmath.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

/**
 * calculate the rotation matrix R.
 */
math::Matrixl math::rmat(math::Box const & box) {
  math::Vecl Rx = box(0) / math::abs(box(0));
  math::Vecl Ry_aux = box(1)
          - math::dot(box(0), box(1)) * box(0)
          / (math::abs(box(0)) * math::abs(box(1)));
  math::Vecl Ry = Ry_aux / math::abs(Ry_aux);
  math::Vecl Rz = math::cross(Rx, Ry);
  math::Matrixl R(Rx, Ry, Rz);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (fabs(R(i, j)) <= math::epsilon)
        R(i, j) = 0.0;
  return R;
}

/**
 * calculate the transformation matrix S.
 */
math::Matrixl math::smat(math::Box const & box, math::boundary_enum const b) {
  if (b == math::rectangular || b == math::truncoct) {
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(0.0, 1.0, 0.0);
    math::Vecl Sz(0.0, 0.0, 1.0);
    math::Matrixl S(Sx, Sy, Sz);
    return S;
  } else {
    long double alpha, beta, gamma;
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
        if (fabs(S(i, j)) < math::epsilon)
          S(i, j) = 0.0;
    return S;
  }
}

/**
 * calculate the inverse of the transformation matrix S.
 */
math::Matrixl math::sinvmat(math::Box const & box, math::boundary_enum const b) {
  if (b == math::rectangular || b == math::truncoct) {
    math::Vecl Sx(1.0, 0.0, 0.0);
    math::Vecl Sy(0.0, 1.0, 0.0);
    math::Vecl Sz(0.0, 0.0, 1.0);
    math::Matrixl Sinv(Sx, Sy, Sz);
    return Sinv;
  } else {
    long double alpha, beta, gamma;
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
