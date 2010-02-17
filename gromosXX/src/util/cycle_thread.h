/**
 * @file cycle_thread.h
 * Thread, that performs a cycle
 */

#ifndef _CYCLE_THREAD_H
#define	_CYCLE_THREAD_H
#include <pthread.h>
#include "thread.h"

namespace util {

  /**
   * @class CycleThread
   * Class for a Thread, that performs a Cycle, everytime
   * do_cycle() is called
   */
  class CycleThread : public util::Thread {
  public:
    /**
     * Consttuctor
     */
    CycleThread();
    /**
     * Where the cycle (while loop) is placed
     */
    virtual void run();

    /**
     * What to do before the cycle in run()
     */
    virtual void init_run() = 0;

    /**
     * What to do during the cycle in run()
     */
    virtual void cycle() = 0;

    /**
     * What to do after the cycle in run()
     */
    virtual void end_run() = 0;

    /**
     * Initiate a cycle
     */
    void do_cycle();

    /**
     * For the first cycle
     */
    void do_first_cycle();

    /**
     * Called by the derrived class to terminate the thread
     */
    void terminate_cycle();


  private:
    //Mutex to control the threads
    pthread_mutex_t mutex_calc;
    pthread_mutex_t mutex_cycle;
    pthread_cond_t cond_calc;
    pthread_cond_t cond_cycle;
    // Needed to avoid spurious wakeup
    bool do_a_calc;
    bool do_a_cycle;

    /**
     * for the while loop in run()
     */
    bool keeprunning;
  };

}





#endif	/* _CYCLE_THREAD_H */

