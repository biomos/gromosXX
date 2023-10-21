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
 * @file cycle_thread.h
 * Thread, that performs a cycle
 */

#ifndef _CYCLE_THREAD_H
#define	_CYCLE_THREAD_H
#include <pthread.h>
#include "thread.h"
//#include "pthread_barrier.h"

#ifdef __APPLE__

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;

int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count);

int pthread_barrier_destroy(pthread_barrier_t *barrier);

int pthread_barrier_wait(pthread_barrier_t *barrier);

#endif

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
  protected:
    /**
     * Barrier for thread initialization
     */
    pthread_barrier_t barrier_init;
    /**
     * Error integer to communicate problems between threads
     */
    int error;
  private:
    /**
     * Barrier at the beginning of the cyle
     */
    pthread_barrier_t barrier_start;
    /**
     * Barrier at the end of the cyle
     */
    pthread_barrier_t barrier_end;

    /**
     * for the while loop in run()
     */
    bool keeprunning;
  };
}

#endif	/* _CYCLE_THREAD_H */
