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
     * reset sub-timer
     */
    void reset(const std::string & name = "total");
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
     * name of the timer
     */
    std::string m_name;
    /**
     * maximum number of threads tracked by timer
     */
    uint m_max_thread_num;

    /**
     * struct sub-timer
     */
    struct subtimer_struct {
      uint master_thread = std::numeric_limits<uint>::max();
      uint highest_used_thread = 0;
      std::vector<double> start_time = start_time = std::vector<double>(1, 0.00);
      std::vector<double> end_time = start_time = std::vector<double>(1, 0.00);
      double total_walltime = 0.0;               //total runtime for subtimer
      double total_cputime = 0.0;       
      bool double_counted_time = false;   //set true if the subtime was once started while another subtimer was still running
    };
    std::map<std::string, subtimer_struct> m_subtimers;
    /**
     * store which subtimer is actually in use (debug only)
     */
    std::string subtimer_in_use = "none";
    /**
     * bool; if detailed report should be printed
     */
    bool m_detailed_report = false;
    /**
     * store the current step of the simulation
     */    
    uint m_step_index;
  };
}

#endif


  
