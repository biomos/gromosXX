/**
 * @file timing.cc
 */

#include <util/stdheader.h>

#include <stddef.h> /* defines NULL */
#include <sys/time.h>

#include "timing.h"

/*
 * writes the current time in seconds to the argument *d.
 * precision (as advertised on the above url): microseconds.
 */
double util::now()     
{
  struct timeval tp;
  gettimeofday(&tp, NULL);

  return ((double)tp.tv_sec+(1.e-6) * tp.tv_usec);
}
