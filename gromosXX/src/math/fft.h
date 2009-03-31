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

