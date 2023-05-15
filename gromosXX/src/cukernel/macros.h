/**
 * @file macros.h
 * defines some macros, e.g. if single or double precision
 */

#ifndef INCLUDED_CUKERNEL_MACROS_H
#define	INCLUDED_CUKERNEL_MACROS_H

#include <iostream>

#define SINGLE_PREC
//#define DOUBLE_PREC

#ifdef SINGLE_PREC
#define REAL float
#define PREC(x) x ## f
#define VECTOR float3
#define MATRIX float9
#define MAKE_VECTOR(x, y, z) make_float3((x), (y), (z))
#define HOST_NEW_POS gpu_stat->host_new_pos
#define HOST_OLD_POS gpu_stat->host_old_pos
#define DEV_NEW_POS gpu_stat->dev_new_pos
#define DEV_OLD_POS gpu_stat->dev_old_pos
#define DEV_FACTOR gpu_stat->dev_factor
#define DEV_CONST_LENGTH2 gpu_stat->dev_const_length2
#define DEV_TOL gpu_stat->dev_tol
#define DEV_MASS gpu_stat->dev_mass
#define FABS(x) fabsf(x)

#elif defined DOUBLE_PREC
#define REAL double
#define PREC(x) x
#define VECTOR double3
#define MATRIX double9
#define MAKE_VECTOR(x, y, z) make_double3((x), (y), (z))
#define HOST_NEW_POS gpu_stat->host_double_new_pos
#define HOST_OLD_POS gpu_stat->host_double_old_pos
#define DEV_NEW_POS gpu_stat->dev_double_new_pos
#define DEV_OLD_POS gpu_stat->dev_double_old_pos
#define DEV_FACTOR gpu_stat->dev_double_factor
#define DEV_CONST_LENGTH2 gpu_stat->dev_double_const_length2
#define DEV_TOL gpu_stat->dev_double_tol
#define DEV_MASS gpu_stat->dev_double_mass
#define FABS(x) fabs(x)
#else
std::cout << "Define precisions!" << std::endl;
#endif

/* How to perform double->float operations on CPU */
//#define FP_CAST float // the old way, float=float*double // just for the consistency check
#define FP_CAST double // consistent, higher precision float=float(double*double) // preferred

#endif	/* _MACROS_H */

