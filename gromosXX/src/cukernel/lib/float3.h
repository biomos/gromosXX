/**
 * @file float3.h 
 * 3D vector operations
 */

#ifndef HOSTDEVICE
#error "Don't include float3.h without defining HOSTDEVICE"
#else

#include <cstdlib>

/**
 * calculates the nearest image for cubic periodic boundary conditions
 * @param[in] i first position
 * @param[in] j second position
 * @param[in] box parameters as float3. Where x is half box, y is the box length, and z is the inverted box length
 * @return the nearest image distance
 */
HOSTDEVICE float3 nearestImage(const float3 & i, const float3 & j, const float3 & box_param_x, const float3 & box_param_y, const float3 & box_param_z) {
  float3 r;

  r.x = i.x - j.x;
  if (fabs(r.x) > box_param_x.x)
    r.x -= box_param_x.y * rintf(r.x * box_param_x.z);

  r.y = i.y - j.y;
  if (fabs(r.y) > box_param_y.x)
    r.y -= box_param_y.y * rintf(r.y * box_param_y.z);

  r.z = i.z - j.z;
  if (fabs(r.z) > box_param_z.x)
    r.z -= box_param_z.y * rintf(r.z * box_param_z.z);

  return r;
}

/**
 * calculates the scalar (dot) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
HOSTDEVICE float dot(const float3 & a, const float3 & b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * calculates the vector (cross) product of two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return the dot product
 */
HOSTDEVICE float3 cross(const float3 & a, const float3 & b) {
  return make_float3(a.y*b.z - a.z*b.y,
                     a.z*b.x - a.x*b.z,
                     a.x*b.y - a.y*b.x);
}

/**
 * calculates the squared length of a vector
 * @param a the vector
 * @return squared length
 */
HOSTDEVICE float abs2(const float3 & a) {
  return dot(a,a);
}

/**
 * calculates the length of a vector
 * @param a the vector
 * @return length
 */
HOSTDEVICE float abs(const float3 & a) {
  return sqrtf(abs2(a));
}

/**
 * operator to subtract two vectors
 * @param[in] a first vector
 * @param[in] b second vector
 * @return a-b
 */
HOSTDEVICE float3 operator-(const float3 & a, const float3 & b) {
  return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

/**
 * operator to add two vectors
 * @param[in] a first vector
 * @param[out] b second vector
 * @return a+b
 */
HOSTDEVICE float3 operator+(const float3 & a, const float3 & b) {
  return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

/**
 * negates a vector
 * @param[in] a the vector
 * @return -a
 */
HOSTDEVICE float3 operator-(const float3 & a) {
  return make_float3(-a.x,-a.y,-a.z);
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
HOSTDEVICE float3 operator*(const float3 & a, float b) {
  return make_float3(a.x*b, a.y*b, a.z*b);
}

/** 
 * scales a vector (multiplication with a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a*b
 */
HOSTDEVICE float3 operator*(float b, const float3 & a) {
  return make_float3(a.x*b, a.y*b, a.z*b);
}

/** 
 * scales a vector (division by a scalar)
 * @param[in] a the vector
 * @param[in] b the scalar
 * @return a/b
 */
HOSTDEVICE float3 operator/(const float3 & a, float b) {
  b = 1.0f / b;
  return a*b;
}

#endif

