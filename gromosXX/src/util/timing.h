/**
 * @file timing.h
 * timing routines.
 */

#ifndef HAVE_INCLUDED_TIMING_H
#define HAVE_INCLUDED_TIMING_H

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
     * Constructor
     */
    Algorithm_Timer() : m_name("anonymous timer") {}
    /**
     * Constructor with name
     */
    Algorithm_Timer(const std::string & name) : m_name(name) {}
    /**
     * start sub timer
     */
    void start(const std::string & name = "total");
    /**
     * stop sub timer
     */
    void stop(const std::string & name = "total");
    /**
     * reset sub timer
     */
    void reset(const std::string & name = "total");
    /**
     * get the total time of a sub timer
     */
    double get_total(const std::string & name = "total");
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
     * the current subtimes
     */
    std::map<std::string, double> m_current_subtimes;
    /** 
     * the total subtimes
     */
    std::map<std::string, double> m_total_subtimes;
  };
}

#endif


  
