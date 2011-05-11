/**
 * @file timing.cc
 */

#include "../stdheader.h"

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

/**
 * starts a (sub)timer
 */
void util::Algorithm_Timer::start(const std::string & name) {
  if (m_current_subtimes.find(name) == m_current_subtimes.end())
    m_total_subtimes[name] = 0.0; // set totals to zero if not exists
  
  m_current_subtimes[name] = util::now();
}
/**
 * stops a (sub)timer
 */
void util::Algorithm_Timer::stop(const std::string & name) {
  if (m_current_subtimes.find(name) == m_current_subtimes.end())
    throw std::runtime_error(name + " timer was stoppped but not stared.");
  
  m_total_subtimes[name] += util::now() - m_current_subtimes[name];
}
/**
 * resets a (sub) timer
 */
void util::Algorithm_Timer::reset(const std::string & name) {
  m_total_subtimes[name] = 0.0;
}
/**
 * gets the total time for a (sub)timer
 */
double util::Algorithm_Timer::get_total(const std::string & name) {
  if (m_total_subtimes.find(name) == m_total_subtimes.end())
    throw std::runtime_error(std::string("Timer ") + m_name + ": " + name + 
            " timer does not exist.");  
  
  return m_total_subtimes[name];
}
/**
 * print the timings to a stream
 */
void util::Algorithm_Timer::print(std::ostream &os) {
  
  if (m_total_subtimes.empty())
    return;
  os.setf(std::ios::fixed, std::ios::floatfield);
  
  // get the totals to calculate 
  const double total = get_total();

  os.precision(6);
  os << std::left << std::setw(10) << " " << std::left << std::setw(30) << m_name 
     << std::setw(14) << std::right << total << std::endl;
  
  std::map<std::string, double>::const_iterator
      it = m_total_subtimes.begin(),
      to = m_total_subtimes.end();
  
  for(; it != to; ++it) {
    os.precision(6);
    if (it->first == "total") continue;
    os << std::right << std::setw(15) << "             - " 
       << std::left << std::setw(25) << it->first
       // totals and percentage
       << std::setw(14) << std::right << it->second;
    os.precision(2);
    os << std::setw(10) << std::right
       << ((total != 0.0) ? (it->second * 100.0 / total) : 0.0) << "%"
       << std::endl;
  }
  
}
