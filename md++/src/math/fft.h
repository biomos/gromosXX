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
 * @file fft.h
 * a wrapper for FFTW
 */

#ifndef _FFT_H
#define	_FFT_H

// make sure we have the complex number type
#include <complex>

#ifdef XXMPI
#include <mpi.h>
#include <fftw3-mpi.h>
#else
#include <fftw3.h>
#endif

#include <config.h>
// these macros will create the FFTW libraries function names.
// it has to be that complicated due to interpolation of the FFTW_PREFIX
// define in config.h
#define GROMOS_FFTW_CONCAT_X(prefix,name) prefix ## name
#define GROMOS_FFTW_CONCAT(prefix,name) GROMOS_FFTW_CONCAT_X(prefix,name)
#define FFTW3(name) GROMOS_FFTW_CONCAT(FFTW_PREFIX,name)
#endif	/* _FFT_H */

