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
 * @file macros.h
 * defines some macros, e.g. if single or double precision
 */

#ifndef _MACROS_H
#define	_MACROS_H

#include <iostream>

#define SINGLE_PREC
//#define DOUBLE_PREC

#ifdef SINGLE_PREC
#define FL_PT_NUM float
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
#define FL_PT_NUM double
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


#endif	/* _MACROS_H */

