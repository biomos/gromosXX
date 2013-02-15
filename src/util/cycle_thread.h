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
     * Called by the derrived class to terminate the thread
     */
    void terminate_cycle();
  private:
    pthread_barrier_t barrier_start;
    pthread_barrier_t barrier_end;

    /**
     * for the while loop in run()
     */
    bool keeprunning;
  };
}

#endif	/* _CYCLE_THREAD_H */
