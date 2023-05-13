/**
 * @file types.h
 * additional types of general use
 */

#ifndef INCLUDED_TYPES_H
#define INCLUDED_TYPES_H

/* *
 * @struct double3 a structure to hold a vector in double precision
 * this is now standardly available and thus commented out!
 */
// struct double3 {
//   /** coordinate x */ double x;
//   /** coordinate y */ double y;
//   /** coordinate z */ double z;
// };

/**
 * @struct float9 a matrix in single precision
 */
struct float9 {
    /** coordinate xx */ float xx;
    /** coordinate xy */ float xy;
    /** coordinate xx */ float xz;
    /** coordinate yx */ float yx;
    /** coordinate yy */ float yy;
    /** coordinate yz */ float yz;
    /** coordinate zx */ float zx;
    /** coordinate zy */ float zy;
    /** coordinate zz */ float zz;
};

/**
 * @struct float9 a matrix in single precision
 */
struct double9 {
    /** coordinate xx */ double xx;
    /** coordinate xy */ double xy;
    /** coordinate xx */ double xz;
    /** coordinate yx */ double yx;
    /** coordinate yy */ double yy;
    /** coordinate yz */ double yz;
    /** coordinate zx */ double zx;
    /** coordinate zy */ double zy;
    /** coordinate zz */ double zz;
};

#endif /* TYPES_H */
