/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file types.h
 * additional types of general use
 */

#ifndef TYPES_H
#define TYPES_H

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
