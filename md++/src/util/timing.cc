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
 * @file timing.cc
 */

#include "../stdheader.h"
#include "../util/debug.h"
#include "../simulation/simulation.h"

#include <stddef.h> /* defines NULL */
#ifdef OMP
  #include <omp.h>
#endif
#ifdef XXMPI
    #include <mpi.h>
#endif

#ifdef WIN32
  #include <time.h>
  #include <sys/types.h>
  #include <sys/timeb.h>
#else
  #include <sys/time.h>
#endif

#include "timing.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE timing

/*
 * writes the current time in seconds to the argument *d.
 * precision: microseconds.
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
 * Constructor
 */
util::Algorithm_Timer::Algorithm_Timer(const std::string & name) {
  DEBUG(5, "Construct main-timer: " + m_name);
  m_name = name;

  /* 
  *  Evaluate number of available threads
  */
  uint tracked_thread_count = 1;
  #ifdef OMP
    #pragma omp parallel
    {
      if (omp_get_thread_num() == 0){   //Main thread
		    tracked_thread_count = omp_get_num_threads();
        tracked_thread_count += 8;     //Extra threads for GPUs
        tracked_thread_count += 2;     //Additional extra threads
      }
    }
  #endif
  m_max_thread_num = tracked_thread_count;
  DEBUG(10, "Maximum tracked threads in main-timer " + m_name + " : " + std::to_string(m_max_thread_num)); 
}

/**
 * initializes a (sub) timer
 */
void util::Algorithm_Timer::initialize(const std::string & name) {
  DEBUG(5, "Initializes sub-timer " + name + " of main-timer " + m_name + ".");

  std::pair<Subtimer_Container::iterator, bool>
    ret = m_subtimers.emplace(name, Subtimer(m_max_thread_num));
  // assert the insertion took place
  assert(ret.second && "re-initialization of sub-timer is not allowed");
  subtimer_in_use.clear();
}

/**
 * starts a main-timer (The main-timer has to be created before the attached sub-timers)
 */
void util::Algorithm_Timer::start(const simulation::Simulation & sim) {
  
  DEBUG(5, "Start main-timer");

  /* 
  *  Store simulation step
  */
  m_step_index = sim.steps();

  /**
   * Create detailed timing output. Set "@print timings_report" as md++ input parameter
   */  
  m_detailed_report = false;
  //If "@print timings_report" is set as md++ parameter
  if (sim.param().print.timings_report) {

    //Write detailed report for timesteps 900 to 904 (if pairlist skipstep = 5)
    const uint first_report_step = (180 * sim.param().pairlist.skip_step);
    const uint last_report_step = (first_report_step + sim.param().pairlist.skip_step - 1); 
    if (m_step_index >= first_report_step && m_step_index <= last_report_step) {
      DEBUG(10, "Write detailed timing report");
      m_detailed_report = true;
    
      /**
       * POSIX Threads (pthreads) used by CUDA implementation do not report a OpenMP thread-number and are always ID'd as thread 0
       */
      if (sim.param().innerloop.method == simulation::sla_cuda){
        DEBUG(10, "In the detailed report, CUDA threads are misidetified as thread #0");
      }

      if (sim.param().constraint.solvent.algorithm == simulation::constr_gpu_shake){
        DEBUG(10, "In the detailed report, CUDA threads are misidetified as thread #0");
      }
    }
  }

  #ifdef XXMPI
    // If MPI is used, only print out report for rank 0
    if (sim.mpiControl().threadID != 0){
      m_detailed_report = false;
    }
  #endif  

  #ifdef OMP
    if (omp_get_thread_num() != 0){
      //io::messages.add("Using a main-timer (see util/timing.cc) in parallel regions is not tested!","Timing",io::message::notice);
      throw std::runtime_error("Using a main-timer (see util/timing.cc) in parallel regions is not tested!");
    }
  #endif

  start_subtimer("total");    //Each main-timer stores its total runtime in a sub-timer called "total"
}

/**
 * starts a sub-timer
 */
void util::Algorithm_Timer::start_subtimer(const std::string & name) {

  uint thread_id = 0; //ID of the thread

  #ifdef OMP
  #pragma omp critical (timing_system)
  #endif
  {
    /**
     * If the subtimer is called first, initialize the variables
     */

      DEBUG(5, "Start sub-timer " + name);
      if (m_subtimers.count(name) == 0){
        {
          initialize(name);
        }
      }
    
    /**
     * Multithreading is handled differently in OpenMP and MPI
     * 
     * In the OpenMP implementation, the main-timer is created once at simulation setup and 
     * called by the master-thread at the start of each simulation step.
     * All the sub-timer for all OpenMP-threads exist within the one object ceated as main-timer by the master-thread.
     * 
     * In the MPI implementation, one main-timer is created per MPI process. 
     * Each MPI process behaves as single-thread application. 
     */
    #ifdef OMP
      thread_id = omp_get_thread_num();
      if (thread_id >= m_max_thread_num) {
        throw std::runtime_error("Maximum number of threads exceeded! Change util/timing.cc to allow more OpenMP threads\n");
      }
      //Find the highest thread number used by one task and save it.
      const uint highest_thread = m_subtimers.at(name).highest_used_thread;
      m_subtimers.at(name).highest_used_thread = std::max(highest_thread, thread_id);
    #endif

    #ifdef XXMPI
      thread_id = 0;    //The MPI Process has always just one thread. If this is changed in the future, the implementation of the timer has to be changed
      #ifdef OMP
        io::messages.add("Usage of both OpenMP and MPI not yet implemented in util/timing.cc. Timings might be wrong!","Timing",io::message::notice);
      #endif
    #endif

    //subtimer.master_thread gets initialized at max(uint); set it to the lowest thread availiable
    //This assumes that the ID of the master thread is constant. Usually it is zero. It can be higher for GPU threads.
    const uint master_thread = m_subtimers.at(name).master_thread;
    m_subtimers.at(name).master_thread = std::min(master_thread, thread_id);

    m_subtimers.at(name).start_time.at(thread_id) = util::now();
  } // End omp critical

  #ifndef NDEBUG  //DEBUG
  //Enable subtimer-lock if
  //  - It is a subtimer
  //  - If it is masterthread
  #ifdef OMP
  #pragma omp critical (timing_system)
  #endif
  {
    if (name != "total" && thread_id == m_subtimers.at(name).master_thread){ //Enable subtimer-lock only if it is a master-thread

      std::string jobname = name + std::to_string(thread_id);
      
      if (subtimer_in_use.empty()){
            subtimer_in_use = jobname;
      } else {
        m_subtimers.at(name).double_counted_time = true;
        DEBUG(5, "Timer " + jobname + " has started while timer " + subtimer_in_use + " is running.");
      }
    }
  } // End omp critical
  #endif
}

/**
 * stops a main-timer
 */
void util::Algorithm_Timer::stop() {
  DEBUG(5, "Stop main-timer");
  stop_subtimer("total");
  
  if (m_detailed_report){
    print_report(std::cout);
  } 

}

/**
 * stops a sub-timer
 */
void util::Algorithm_Timer::stop_subtimer(const std::string & name) {

  DEBUG(5, "Stop sub-timer");
  uint thread_id = 0;
  #ifdef OMP
      thread_id = omp_get_thread_num();
  #endif

  #ifdef XXMPI
    thread_id = 0;    //The MPI Process has always just one thread. If this is changed in the future, the implementation of the timer has to be changed
    #ifdef OMP
      io::messages.add("Usage of both OpenMP and MPI not yet implemented in util/timing.cc. Timings might be wrong!","Timing",io::message::notice);
    #endif
  #endif

  #ifdef OMP
  #pragma omp critical (timing_system)
  #endif
  {  /*
    *  Check if the timer has been stopped without beeing in the subtimer map.
    *  This should never happen!!
    *  It would mean, that the subtimer was not yet initialized and we have a serious race condition
    *
    * calling .count() without lock is thread-unsafe
    * existence of 'name' means initialized subtimer, however
    * adding another subtimer simultaneously by another thread
    * may lead to rehashing, that temporarily invalidates
    * the .count() call
    */ 
    if (m_subtimers.count(name) == 0) {  
      {
        //io::messages.add("Sub-timer " + name + " was stopped that has not been initialized before!","Timing",io::message::notice);
        //DEBUG(5, "Sub-timer " + name + " was stopped that has not been initialized before!");
        throw std::runtime_error("Sub-timer " + name + " was stopped that has not been initialized before! (Simulation step: " + std::to_string(m_step_index) + ")\n");
      }
    }
    m_subtimers.at(name).end_time.at(thread_id) = util::now();

    if (thread_id == m_subtimers.at(name).master_thread){
      m_subtimers.at(name).total_walltime += (m_subtimers.at(name).end_time.at(thread_id) - m_subtimers.at(name).start_time.at(thread_id));
    }
    m_subtimers.at(name).total_cputime += (m_subtimers.at(name).end_time.at(thread_id) - m_subtimers.at(name).start_time.at(thread_id));

  /*
  * Disable subtimer_in_use flag after stopping subtimer
  */
  #ifndef NDEBUG
    std::string jobname = name + std::to_string(thread_id);
    if (name != "total" && subtimer_in_use == jobname){
      subtimer_in_use.clear();
    }    
  #endif
  }
}


/**
 * gets the total time for a (sub)timer
 */
double util::Algorithm_Timer::get_total(const std::string & name) {

  if (m_subtimers.count(name) == 0){
    throw std::runtime_error(std::string("Timer ") + m_name + ": " + name + 
            " timer does not exist.");  
  }
  
  //Sanity check: Throw error if runtime for a subtimer is larger than 2 years
  /*if (std::abs(m_subtimers.at(name).total_walltime) > 63072000.0){
    //io::messages.add("Total runtime of timer " + name + " failed sanity check!","Timing",io::message::notice);
    //DEBUG(5, "Total runtime of timer " + name + " failed sanity check!");
    throw std::runtime_error("Total runtime of timer " + name + " failed sanity check!"); //DELETE BEFORE MERGE
  }*/

  return m_subtimers.at(name).total_walltime;
}
 
/**
 * print detailed report to file
 */
void util::Algorithm_Timer::print_report(std::ostream & os) {

  os << "DETAILED_TIMING" << std::endl;
  for(auto& subtimer : m_subtimers) {
    //Print start time
    os << m_name;
    os << ";" << subtimer.first;
    os << ";start";
    os << std::fixed << std::setprecision(6);
    for(uint i = 0; i<=subtimer.second.highest_used_thread; i++){
      os << "-" << subtimer.second.start_time.at(i);
    }    
    os << ";end";
    os << std::fixed << std::setprecision(6);
    for(uint i = 0; i<=subtimer.second.highest_used_thread; i++){
      os << "-" << subtimer.second.end_time.at(i);
    }
    os << std::endl;

    //Delete data after writing
    #ifdef OMP
    #pragma omp critical (timing_system)
    #endif
    {
      subtimer.second.start_time = std::vector<double>(m_max_thread_num, 0.00);
      subtimer.second.end_time = std::vector<double>(m_max_thread_num, 0.00);
      subtimer.second.highest_used_thread = 0;
    }
  }
  os << "END" << std::endl;
}

/**
 * print the timings to a stream
 */
void util::Algorithm_Timer::print(std::ostream &os) {
  
  if (m_subtimers.empty())
    return;
  os.setf(std::ios::fixed, std::ios::floatfield);
  
  // get the totals to calculate 
  const double total = get_total();

  os.precision(3);
  os << std::left << std::setw(10) << " " << std::left << std::setw(30) << m_name 
     << std::setw(14) << std::right << total;
     
  if (m_subtimers.size() > 1){
    os << std::setw(21) << "100.00%" << std::endl;
  } else {
    os << std::endl;
  }
  
  Subtimer_Container::const_iterator
      it = m_subtimers.begin(),
      to = m_subtimers.end();
  
  // Print the subtimers
  double accounted_time = 0.0;
  for(; it != to; ++it) {
    os.precision(3);
    if (it->first == "total") continue;
    os << std::right << std::setw(15) << "             - " 
       << std::left << std::setw(35) << it->first
       // totals and percentage
       << std::setw(14) << std::right << it->second.total_walltime;
    os.precision(2);
    os << std::setw(10) << std::right
       << ((total != 0.0) ? (it->second.total_walltime * 100.0 / total) : 0.0) << "%";
    if (it->second.double_counted_time) os << "      [includes double counted time]"; 
    os << std::endl;
    accounted_time += it->second.total_walltime;
  }
  
  if (accounted_time != 0.0){ //Print remaining time to 100%
    double unaccounted_time = total-accounted_time;
    os.precision(3);
    os << std::right << std::setw(15) << "             - " 
       << std::left << std::setw(35) << "unaccounted time"
       << std::setw(14) << std::right << (unaccounted_time);
    os.precision(2);
    os << std::setw(10) << std::right
       << ((total != 0.0) ? (100.0 * (unaccounted_time) / total) : 0.0) << "%";
    if (unaccounted_time < 0.0) os << "      [some timings above may include double counted time]"; 
    os << std::endl;
  }
  
}
