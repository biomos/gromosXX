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
 * @file timing.h
 * timing routines.
 */

#ifndef HAVE_INCLUDED_TIMING_H
#define HAVE_INCLUDED_TIMING_H

namespace simulation{
	class Simulation;
}

namespace util
{
  /**
   * return the time now.
   */
  double now();
  
  /**
   * @class Algorithm_Timer
   * a class to do timings of algorithms
   */
  class Algorithm_Timer {
  public:
    /**
     * Constructor with name
     */
    Algorithm_Timer(const std::string & name);
    /**
     * start main-timer
     */
    void start(const simulation::Simulation & sim);
    /**
     * start sub-timer
     */
    void start_subtimer(const std::string & name = "total");
    /**
     * stop main-timer
     */
    void stop();
    /**
     * stop sub-timer
     */
    void stop_subtimer(const std::string & name = "total");
    /**
     * initialize sub-timer
     */
    void initialize(const std::string & name = "total");
    /**
     * get the total time of a sub-timer
     */
    double get_total(const std::string & name = "total");
    /**
    * print detailed report to omd-file
    */
    void print_report(std::ostream & out = std::cout);
    /**
     * print timing results to stream
     */
    void print(std::ostream & out = std::cout);
    /**
     * accessor to name
     */
    void name(const std::string & name) {
      m_name = name;
    }
    /**
     * accessor to name
     */
    std::string & name() {
      return m_name;
    }
    /**
     * const accessor to name
     */
    const std::string & name() const {
      return m_name;
    }    
    
  
  protected:
    /**
     * name of the maintimer timer
     */
    std::string m_name;
    /**
     * maximum number of threads tracked by timer
     */
    uint m_max_thread_num;
    /**
     * store the current step of the simulation
     */    
    uint m_step_index;
    /**
     * bool; if detailed report should be printed
     */
    bool m_detailed_report = false;

    /**
     * sub-timer class
     */
    class Subtimer {
      public:
        /**
         * vector of start times - one per thread
         */
        std::vector<double> start_time;
        /**
         * vector of end times - one per thread
         */
        std::vector<double> end_time;
        /**
         * total walltime for subtimer
         */
        double total_walltime = 0.0;
        /**
         * total cputime for subtimer
         */
        double total_cputime = 0.0;
        /**
         * ID of the master thread
         */  
        uint master_thread = std::numeric_limits<uint>::max();
        /**
         * highest ID of the used thread
         */
        uint highest_used_thread = 0;
        /**
         * true if the subtimer was started while another subtimer was still running
         */
        bool double_counted_time = false;

        /**
         * Constructor
         */
        Subtimer(const uint thread_count = 1){
          start_time = std::vector<double>(thread_count, 0.0);
          end_time = std::vector<double>(thread_count, 0.0);
        }
    };
    typedef std::unordered_map<std::string, Subtimer> Subtimer_Container;
    /**
     * a map of subtimers, one subtimer per algorithm name
     */
    Subtimer_Container m_subtimers;
    /**
     * store which subtimer is actually in use (debug only)
     */
    std::string subtimer_in_use;
  };
}

#endif


  
