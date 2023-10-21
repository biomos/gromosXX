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
 * @file cycle_thread.cc
 * Thread performin a cycle
 */

#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include "../stdheader.h"
#include "../util/debug.h"
//#include "../util/pthread_barrier.h"
#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

#include "cycle_thread.h"

#ifdef __APPLE__
int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count)
{
    if(count == 0)
    {
        errno = EINVAL;
        return -1;
    }
    if(pthread_mutex_init(&barrier->mutex, 0) < 0)
    {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }
    barrier->tripCount = count;
    barrier->count = 0;

    return 0;
}

int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }
    else
    {
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

#endif

/**
 * Constructor
 */
util::CycleThread::CycleThread()
  : error(0) {
  // For the while loop in run()
  keeprunning = true;
  // Initialize all the barriers
  pthread_barrier_init(&barrier_init, NULL, 2);
  pthread_barrier_init(&barrier_start, NULL, 2);
  pthread_barrier_init(&barrier_end, NULL, 2);
  // The next cycle will be the first
  DEBUG(15, "CycleThread: Initialized CycleThread variables")
}

/**
 * Called by the derrived class to terminate the thread
 */
void util::CycleThread::terminate_cycle() {
  // Terminate while loop in run()
  keeprunning = false;
  // Make sure it terminates by doing one more cycle
  do_cycle();
  // wait for thread to finish
  wait();
  pthread_barrier_destroy(&barrier_init);
  pthread_barrier_destroy(&barrier_start);
  pthread_barrier_destroy(&barrier_end);
  DEBUG(15, "CycleThread: Destroyed all CycleThread variables")
}

/**
 * Where the cycle (while loop) is placed
 */
void util::CycleThread::run() {
  init_run();
    pthread_barrier_wait(&barrier_init);
  while (keeprunning) {
    // wait till told to start
    pthread_barrier_wait(&barrier_start);
    cycle();
    // wait till finished
    pthread_barrier_wait(&barrier_end);
  }
  end_run();
}

/**
 * Initiate a cycle
 */
void util::CycleThread::do_cycle() {
  // start the cycle 
  pthread_barrier_wait(&barrier_start);
  // wait for the cycle to finish
  pthread_barrier_wait(&barrier_end);
}
