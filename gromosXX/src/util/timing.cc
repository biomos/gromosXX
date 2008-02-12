/**
 * @file timing.cc
 */

#include <stdheader.h>

#include <stddef.h> /* defines NULL */
#ifdef WIN32

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#else

#include <sys/time.h>

#endif

#include "timing.h"

/*
 * writes the current time in seconds to the argument *d.
 * precision (as advertised on the above url): microseconds.
 */
double util::now()     
{
#ifndef WIN32
	struct timeval tp;
  gettimeofday(&tp, NULL);

  return ((double)tp.tv_sec+(1.e-6) * tp.tv_usec);
#else
	struct __timeb64 tstruct;
	_ftime64( &tstruct );
	return double(tstruct.time + (1.e-6) * tstruct.millitm);
#endif
}
